/**
 * ARCANA API bridge — maps frontend calls to FastAPI backend endpoints.
 * Loaded by index.html before any dashboard code runs.
 */

const ARCANA_API = (() => {
  // Auto-detect the backend origin (same host, port 8765 if served differently)
  const BASE = window.location.port === '8765'
    ? window.location.origin
    : `${window.location.protocol}//${window.location.hostname}:8765`;

  async function _get(path, params = {}) {
    const url = new URL(BASE + path);
    Object.entries(params).forEach(([k, v]) => { if (v !== undefined && v !== null && v !== '') url.searchParams.set(k, v); });
    const r = await fetch(url.toString());
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

  return {
    BASE,

    /** Returns true if the backend is reachable */
    async ping() {
      try {
        await fetch(BASE + '/api/config', { signal: AbortSignal.timeout(3000) });
        return true;
      } catch { return false; }
    },

    /** Start a new analysis run */
    async startRun(params = {}) {
      return _post('/api/run', params);
    },

    /** List past runs */
    async listRuns() {
      return _get('/api/runs');
    },

    /** Get a single run */
    async getRun(runId) {
      return _get(`/api/run/${runId}`);
    },

    /**
     * Open a WebSocket for live run progress.
     * handlers: { onMessage, onDone, onError }
     */
    connectRunWS(runId, handlers = {}) {
      const wsBase = BASE.replace(/^http/, 'ws');
      const ws = new WebSocket(`${wsBase}/ws/run/${runId}`);
      ws.onmessage = e => handlers.onMessage?.(JSON.parse(e.data));
      ws.onclose   = ()  => handlers.onDone?.();
      ws.onerror   = e  => handlers.onError?.(e);
      return ws;
    },

    /** Fetch DE results for a comparison */
    async getDEResults(comparison = '', lncrnaOnly = false) {
      return _get('/api/results/de', { comparison, lncrna_only: lncrnaOnly || undefined });
    },

    /** List available DE comparisons */
    async listDEComparisons() {
      return _get('/api/results/de/comparisons');
    },

    /** Fetch hub genes */
    async getHubGenes({ lncrnaOnly = false, module = '', limit = 50 } = {}) {
      return _get('/api/results/hub_genes', { lncrna_only: lncrnaOnly || undefined, module: module || undefined, limit });
    },

    /** Fetch focus modules */
    async getFocusModules() {
      return _get('/api/results/focus_modules');
    },

    /** Fetch module–trait correlation data */
    async getModuleTrait() {
      return _get('/api/results/module_trait');
    },

    /** Fetch enrichment results */
    async getEnrichment(comparison = '') {
      return _get('/api/results/enrichment', { comparison: comparison || undefined });
    },

    /** Fetch the run summary (module/gene counts etc.) */
    async getSummary() {
      return _get('/api/results/summary');
    },

    /** List available figure files */
    async listFigures() {
      return _get('/api/figures');
    },

    /** Upload a counts matrix file */
    async uploadCounts(file) {
      const fd = new FormData();
      fd.append('file', file);
      const r = await fetch(BASE + '/api/upload/counts', { method: 'POST', body: fd });
      if (!r.ok) throw new Error(`Upload counts → ${r.status}`);
      return r.json();
    },

    /** Upload a metadata file */
    async uploadMetadata(file) {
      const fd = new FormData();
      fd.append('file', file);
      const r = await fetch(BASE + '/api/upload/metadata', { method: 'POST', body: fd });
      if (!r.ok) throw new Error(`Upload metadata → ${r.status}`);
      return r.json();
    },

    /** Search indexed result documents */
    async searchDocs(query) {
      return _get('/api/docs/search', { q: query });
    },

    /** Check AI (Ollama) status */
    async getAIStatus() {
      return _get('/api/ai/status');
    },

    /** Send a message to ARCANA-AI */
    async chatAI(message, history = []) {
      return _post('/api/ai/chat', { message, history });
    },
  };
})();
