'use strict';
// Preload runs in the renderer with Node access bridged via contextBridge.
// Currently unused — kept as an extension point.
const { contextBridge } = require('electron');

contextBridge.exposeInMainWorld('arcanaDesktop', {
  version: process.versions.electron,
});
