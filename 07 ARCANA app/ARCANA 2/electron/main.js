'use strict';

const { app, BrowserWindow, shell, dialog, Menu } = require('electron');
const path  = require('path');
const http  = require('http');
const { spawn } = require('child_process');
const fs    = require('fs');

// ── Configuration ─────────────────────────────────────────────────────────────

const PORT     = 8765;
const BASE_URL = `http://localhost:${PORT}`;

// When packaged, the ARCANA project lives in Resources/arcana/
// When running in dev, it's the parent directory of electron/
const IS_PACKAGED = app.isPackaged;
const ARCANA_ROOT = IS_PACKAGED
  ? path.join(process.resourcesPath, 'arcana')
  : path.resolve(__dirname, '..');

// Python venv uvicorn path
const VENV_PYTHON = IS_PACKAGED
  ? path.join(ARCANA_ROOT, 'venv', 'bin', 'uvicorn')
  : path.join(ARCANA_ROOT, 'venv', 'bin', 'uvicorn');

// ── State ─────────────────────────────────────────────────────────────────────

let backendProcess = null;
let mainWindow     = null;
let splashWindow   = null;

// ── Start FastAPI backend ─────────────────────────────────────────────────────

function startBackend() {
  if (!fs.existsSync(VENV_PYTHON)) {
    dialog.showErrorBox(
      'Setup required',
      `ARCANA venv not found at:\n${VENV_PYTHON}\n\nRun install.sh first.`
    );
    app.quit();
    return;
  }

  console.log('[ARCANA] Starting backend…');
  backendProcess = spawn(VENV_PYTHON, [
    'app.backend.main:app',
    '--host', '127.0.0.1',
    '--port', String(PORT),
  ], {
    cwd: ARCANA_ROOT,
    env: { ...process.env, PYTHONUNBUFFERED: '1' },
  });

  backendProcess.stdout.on('data', d => process.stdout.write(`[backend] ${d}`));
  backendProcess.stderr.on('data', d => process.stderr.write(`[backend] ${d}`));
  backendProcess.on('exit', code => {
    console.log(`[ARCANA] Backend exited (${code})`);
  });
}

// ── Poll until the backend is ready ──────────────────────────────────────────

function waitForBackend(maxWaitMs = 30000, intervalMs = 400) {
  return new Promise((resolve, reject) => {
    const deadline = Date.now() + maxWaitMs;
    const check = () => {
      http.get(`${BASE_URL}/api/config`, res => {
        if (res.statusCode === 200) { resolve(); return; }
        retry();
      }).on('error', retry);
    };
    const retry = () => {
      if (Date.now() > deadline) { reject(new Error('Backend startup timeout')); return; }
      setTimeout(check, intervalMs);
    };
    check();
  });
}

// ── Create splash window ──────────────────────────────────────────────────────

function createSplash() {
  splashWindow = new BrowserWindow({
    width: 480,
    height: 280,
    frame: false,
    transparent: true,
    resizable: false,
    alwaysOnTop: true,
    skipTaskbar: true,
    webPreferences: { contextIsolation: true },
  });
  splashWindow.loadFile(path.join(__dirname, 'splash.html'));
  splashWindow.center();
}

// ── Create main window ────────────────────────────────────────────────────────

function createMainWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 900,
    minHeight: 600,
    show: false,
    titleBarStyle: 'hiddenInset',   // Mac: traffic lights overlap title bar
    backgroundColor: '#0a0a0f',
    icon: path.join(__dirname, 'assets', 'icon.png'),
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
    },
  });

  // Build the application menu
  const menu = Menu.buildFromTemplate([
    {
      label: 'ARCANA',
      submenu: [
        { role: 'about' },
        { type: 'separator' },
        {
          label: 'Open Config File',
          accelerator: 'CmdOrCtrl+,',
          click: () => shell.openPath(path.join(ARCANA_ROOT, 'config', 'project.yaml')),
        },
        {
          label: 'Open Results Folder',
          click: () => shell.openPath(path.join(ARCANA_ROOT, 'results')),
        },
        { type: 'separator' },
        { role: 'quit' },
      ],
    },
    {
      label: 'View',
      submenu: [
        { role: 'reload' },
        { role: 'forceReload' },
        { role: 'toggleDevTools' },
        { type: 'separator' },
        { role: 'resetZoom' },
        { role: 'zoomIn' },
        { role: 'zoomOut' },
        { type: 'separator' },
        { role: 'togglefullscreen' },
      ],
    },
    {
      label: 'Edit',
      submenu: [
        { role: 'undo' }, { role: 'redo' },
        { type: 'separator' },
        { role: 'cut' }, { role: 'copy' }, { role: 'paste' },
        { role: 'selectAll' },
      ],
    },
  ]);
  Menu.setApplicationMenu(menu);

  mainWindow.loadURL(BASE_URL);

  mainWindow.once('ready-to-show', () => {
    if (splashWindow && !splashWindow.isDestroyed()) {
      splashWindow.destroy();
      splashWindow = null;
    }
    mainWindow.show();
  });

  mainWindow.on('closed', () => { mainWindow = null; });
}

// ── App lifecycle ─────────────────────────────────────────────────────────────

app.whenReady().then(async () => {
  createSplash();
  startBackend();

  try {
    await waitForBackend();
    createMainWindow();
  } catch (err) {
    if (splashWindow && !splashWindow.isDestroyed()) splashWindow.destroy();
    dialog.showErrorBox('Startup failed', `Backend did not start in time.\n\n${err.message}`);
    app.quit();
  }
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});

app.on('activate', () => {
  if (mainWindow === null) createMainWindow();
});

app.on('will-quit', () => {
  if (backendProcess) {
    console.log('[ARCANA] Shutting down backend…');
    backendProcess.kill('SIGTERM');
  }
});
