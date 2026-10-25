#!/usr/bin/env python3
"""
ARCANA Desktop Launcher
Opens ARCANA in a native Mac window (no browser needed).
Run: python launch_app.py
"""

import subprocess
import sys
import time
import threading
import urllib.request
from pathlib import Path

ROOT   = Path(__file__).resolve().parent
VENV   = ROOT / "venv" / "bin"
PORT   = 8765
URL    = f"http://localhost:{PORT}"

def start_backend():
    proc = subprocess.Popen(
        [str(VENV / "uvicorn"), "app.backend.main:app",
         "--host", "127.0.0.1", "--port", str(PORT)],
        cwd=str(ROOT),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return proc

def wait_for_backend(timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            urllib.request.urlopen(f"{URL}/api/config", timeout=1)
            return True
        except Exception:
            time.sleep(0.4)
    return False

def main():
    import webview

    # Start the FastAPI backend in background
    backend = start_backend()

    # Show a loading window while backend boots
    splash = webview.create_window(
        "ARCANA — Starting…",
        html="""
        <style>
          body { margin:0; background:#0d0d14; display:flex; align-items:center;
                 justify-content:center; height:100vh;
                 font-family:-apple-system,sans-serif; color:#e8d5a3; }
          .box { text-align:center; }
          h2 { font-size:1.6rem; margin-bottom:0.4rem; }
          h2 em { color:#c97d7d; font-style:italic; font-weight:300; }
          p { color:rgba(255,255,255,0.4); font-size:0.8rem;
              letter-spacing:0.08em; text-transform:uppercase; margin-bottom:2rem; }
          .bar { width:240px; height:3px; background:rgba(255,255,255,0.06);
                 border-radius:2px; overflow:hidden; margin:0 auto; }
          .fill { height:100%; background:linear-gradient(90deg,#c97d7d,#e8d5a3);
                  animation:go 1.8s ease-in-out infinite; transform-origin:left; }
          @keyframes go { 0%{transform:scaleX(0)} 70%{transform:scaleX(.9)} 100%{transform:scaleX(1);opacity:0} }
        </style>
        <div class="box">
          <h2>ARCANA / <em>lncRNA</em></h2>
          <p>Co-expression Analysis Platform</p>
          <div class="bar"><div class="fill"></div></div>
        </div>""",
        width=460, height=260, frameless=True, on_top=True,
    )

    def open_main():
        ready = wait_for_backend()
        if not ready:
            webview.windows[0].title = "ARCANA — Backend failed to start"
            return
        # Destroy splash and open main window
        splash.destroy()
        win = webview.create_window(
            "ARCANA",
            url=URL,
            width=1440, height=900,
            min_size=(900, 600),
        )

    threading.Thread(target=open_main, daemon=True).start()

    webview.start(debug=False)

    # Clean up backend when window closes
    backend.terminate()

if __name__ == "__main__":
    main()
