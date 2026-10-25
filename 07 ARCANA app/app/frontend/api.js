/**
 * ARCANA API client
 * Connects the frontend UI to the FastAPI backend on localhost:8000
 * Falls back gracefully when the backend is not running (demo mode).
 */

const ARCANA_API = (() => {
  const BASE = 'http://localhost:8000';
  let _ws = null;
  let _wsCallbacks = {};

  async function _get(path) {
    const r = await fetch(BASE + path);
    if (!r.ok) throw new Error(`GET ${path} → ${r.status}`);
    return r.json();
  }

  async function _post(path, body) {
    const r = await fetch(BASE + path, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    if (!r.ok) throw new Error(`POST ${path} → ${r.status}`);
    return r.json();
  }

  // ── Pipeline ──────────────────────────────────────────────────────────────

  async function startRun({ name = 'New run', steps = null, params = {} } = {}) {
    const body = { name, params };
    if (steps) body.steps = steps;
    return _post('/api/run', body);
  }

  function connectRunWS(runId, { onLog, onStepStart, onStepDone, onError, onComplete } = {}) {
    if (_ws) _ws.close();
    _ws = new WebSocket(`ws://localhost:8000/ws/run/${runId}`);
    _ws.onmessage = (e) => {
      const msg = JSON.parse(e.data);
      switch (msg.type) {
        case 'log':        onLog?.(msg); break;
        case 'step_start': onStepStart?.(msg); break;
        case 'step_done':  onStepDone?.(msg); break;
        case 'step_error': onError?.(msg); break;
        case 'run_complete': onComplete?.(msg); break;
        case 'run_failed':   onError?.(msg); break;
      }
    };
    _ws.onerror = (e) => onError?.({ type: 'ws_error', error: e });
    return _ws;
  }

  async function listRuns()   { return _get('/api/runs'); }
  async function getRun(id)   { return _get(`/api/run/${id}`); }
  async function listSteps()  { return _get('/api/steps'); }

  // ── Results ───────────────────────────────────────────────────────────────

  async function getDEResults(comparison = 'T2D_vs_ND', lncrnaOnly = false) {
    return _get(`/api/results/de?comparison=${comparison}&lncrna_only=${lncrnaOnly}`);
  }
  async function getHubGenes({ lncrnaOnly = false, module = '', limit = 200 } = {}) {
    return _get(`/api/results/hub_genes?lncrna_only=${lncrnaOnly}&module=${module}&limit=${limit}`);
  }
  async function getFocusModules()  { return _get('/api/results/focus_modules'); }
  async function getModules()       { return _get('/api/results/modules'); }
  async function getModuleTrait()   { return _get('/api/results/module_trait'); }
  async function getMLMetrics()     { return _get('/api/results/ml'); }
  async function getEnrichment()    { return _get('/api/results/enrichment'); }
  async function listFigures()      { return _get('/api/figures'); }
  async function figureURL(path)    { return `${BASE}/api/figure?path=${encodeURIComponent(path)}`; }

  // ── File upload ───────────────────────────────────────────────────────────

  async function uploadCounts(file) {
    const fd = new FormData();
    fd.append('file', file);
    const r = await fetch(`${BASE}/api/upload/counts`, { method: 'POST', body: fd });
    return r.json();
  }

  async function uploadMetadata(file) {
    const fd = new FormData();
    fd.append('file', file);
    const r = await fetch(`${BASE}/api/upload/metadata`, { method: 'POST', body: fd });
    return r.json();
  }

  // ── Document search ───────────────────────────────────────────────────────

  async function searchDocs(query)  { return _get(`/api/docs/search?q=${encodeURIComponent(query)}`); }
  async function indexDocs(folder)  { return _post('/api/docs/index', {}); }

  // ── Health check ──────────────────────────────────────────────────────────

  async function ping() {
    try {
      const r = await fetch(`${BASE}/api/steps`, { signal: AbortSignal.timeout(1000) });
      return r.ok;
    } catch { return false; }
  }

  return {
    startRun, connectRunWS, listRuns, getRun, listSteps,
    getDEResults, getHubGenes, getFocusModules, getModules, getModuleTrait, getMLMetrics, getEnrichment,
    listFigures, figureURL, uploadCounts, uploadMetadata,
    searchDocs, indexDocs, ping,
    BASE,
  };
})();

// ── Auto-connect indicator ─────────────────────────────────────────────────────

(async () => {
  const alive = await ARCANA_API.ping();
  const badge = document.querySelector('.nav-badge');
  if (badge) {
    if (alive) {
      badge.textContent = 'Pipeline v2.1 · connected';
      badge.style.color = 'var(--sage, #7aab8a)';
    } else {
      badge.textContent = 'Backend offline · demo mode';
      badge.style.color = 'var(--muted, #8892aa)';
    }
  }

  // Hydrate run history from backend if available
  if (alive) {
    try {
      const runs = await ARCANA_API.listRuns();
      window.__ARCANA_RUNS__ = runs;
    } catch {}
  }
})();
