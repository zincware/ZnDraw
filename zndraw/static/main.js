import { Cache } from './pycom/Cache.js';
import { World } from './World/World.js';
import { setUIEvents } from './UI/UI.js';
import { initJSONEditor } from './UI/json_editor.js';

function setupSocket() {
  const socket = io();
  socket.on('connect', () => {
    socket.emit('connection', { data: "I'm connected!" });
  });
  return socket;
}

function setupUpload(socket) {
  const upload_form = document.getElementById('uploadForm');
  const upload_input = document.getElementById('fileInput');
  const submit_button = document.getElementById('uploadBtn');

  const file = {
    dom: document.getElementById('fileInput'),
    binary: null,
  };

  const reader = new FileReader();

  // Because FileReader is asynchronous, store its
  // result when it finishes reading the file
  reader.addEventListener('load', () => {
    file.binary = reader.result;
    socket.emit('upload', { content: reader.result, filename: file.dom.files[0].name });
  });

  // At page load, if a file is already selected, read it.
  if (file.dom.files[0]) {
    reader.readAsBinaryString(file.dom.files[0]);
  }

  // If not, read the file once the user selects it.
  file.dom.addEventListener('change', () => {
    if (reader.readyState === FileReader.LOADING) {
      reader.abort();
    }

    reader.readAsBinaryString(file.dom.files[0]);
  });
}

function setupSlider(socket, world) { }

function main() {
  const socket = setupSocket();
  const cache = new Cache(socket);
  const container = document.querySelector('#scene-container');
  const world = new World(container, cache, socket);

  setUIEvents(socket, world);
  initJSONEditor(socket, cache, world);
  world.start();

  // socket dicsconnect event
  socket.on('disconnect', () => {
    window.close();
  });

  document.getElementById('ExitBtn').addEventListener('click', () => {
    fetch('/exit', { method: 'GET' });
  });

  document.getElementById('downloadBtn').addEventListener('click', () => {
    socket.emit('download', { atoms_list: cache.getAllAtoms() }, (data) => {
      const blob = new Blob([data], { type: 'text/csv' });
      const elem = window.document.createElement('a');
      elem.href = window.URL.createObjectURL(blob);
      elem.download = 'trajectory.xyz';
      document.body.appendChild(elem);
      elem.click();
      document.body.removeChild(elem);
    });
  });

  document.getElementById('downloadSelectedBtn').addEventListener('click', () => {
    socket.emit('download', { atoms_list: cache.getAllAtoms(), selection: world.selection.selection }, (data) => {
      const blob = new Blob([data], { type: 'text/csv' });
      const elem = window.document.createElement('a');
      elem.href = window.URL.createObjectURL(blob);
      elem.download = 'trajectory.xyz';
      document.body.appendChild(elem);
      elem.click();
      document.body.removeChild(elem);
    });
  });

  const helpBtn = document.getElementById('HelpBtn');

  helpBtn.addEventListener('mouseover', () => {
    new bootstrap.Collapse(document.getElementById('helpBoxCollapse'), { toggle: false }).show();
  });
  helpBtn.addEventListener('mouseout', () => {
    new bootstrap.Collapse(document.getElementById('helpBoxCollapse'), { toggle: false }).hide();
  });

  socket.emit('atoms:request', null, () => {
    world.setStep(0);
    // disable loading screen
    document.getElementById('atom-spinner').style.display = 'none';
  });

  setupUpload(socket);
}

main();
