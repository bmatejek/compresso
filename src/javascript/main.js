var F = null;

window.onload = function() {
  // init drag and drop
  document.getElementById("body").addEventListener("dragover", onDragOver, false);
  document.getElementById("body").addEventListener("dragleave", onDragLeave, false);
  document.getElementById("body").addEventListener("drop", onDrop, false);
};

function onDragOver(e) {
  // prevent the default handler for traversing drag area
  e.preventDefault();
};

function onDragLeave(e) {
  // prevent the default handler for leaving drag area
  e.preventDefault();
};

function onDrop(e) {
  // don't propagate and prevent default
  e.stopPropagation();
  e.preventDefault();

  // get the dropped files
  var files = e.dataTransfer.files;

  // if anything is wrong with the dropped files, exit
  if (typeof files == "undefined" || files.length == 0) {
    return;
  }

  // create a new compresso file
  F = new SEGMENTATION.File(files[0]);
};
