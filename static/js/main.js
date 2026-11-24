// main.js - small bits of interactivity. The website also works without it.

// Tell the CSS that JavaScript is running (used to hide/show cutoff options)
document.documentElement.classList.add("js");


// ---- 1. Light / dark theme button ----
var themeButton = document.getElementById("themeToggle");
if (themeButton) {
  themeButton.addEventListener("click", function () {
    var current = document.documentElement.getAttribute("data-theme");
    if (!current) {
      current = window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
    }
    var next = current === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    try { localStorage.setItem("theme", next); } catch (e) {}
  });
}


// ---- 2. "Upload file" / "Paste text" tabs ----
var tabs = document.querySelectorAll(".tab");
tabs.forEach(function (tab) {
  tab.addEventListener("click", function () {
    tabs.forEach(function (t) { t.classList.remove("active"); });
    document.querySelectorAll(".tab-panel").forEach(function (p) { p.classList.remove("active"); });
    tab.classList.add("active");
    document.getElementById(tab.dataset.tab).classList.add("active");
  });
});

// If the user had pasted text before an error, open the "Paste text" tab
var pasteBox = document.querySelector("textarea[name='expression_text']");
if (pasteBox && pasteBox.value.trim() !== "") {
  document.querySelector(".tab[data-tab='paste-expr']").click();
}


// ---- 3. Drag-and-drop file box ----
var dropzone = document.querySelector(".dropzone");
if (dropzone) {
  var fileInput = dropzone.querySelector("input[type='file']");
  var fileLabel = dropzone.querySelector(".dz-file");

  fileInput.addEventListener("change", function () {
    if (fileInput.files.length > 0) {
      var file = fileInput.files[0];
      var sizeKb = Math.round(file.size / 1024);
      fileLabel.textContent = "✓ " + file.name + " (" + sizeKb + " KB)";
      dropzone.classList.add("has-file");
    }
  });
  ["dragenter", "dragover"].forEach(function (name) {
    dropzone.addEventListener(name, function () { dropzone.classList.add("dragover"); });
  });
  ["dragleave", "drop"].forEach(function (name) {
    dropzone.addEventListener(name, function () { dropzone.classList.remove("dragover"); });
  });
}


// ---- 4. Show the extra box for the chosen cutoff option ----
function updateCutoffBoxes() {
  var checked = document.querySelector("input[name='cutoff_mode']:checked");
  document.querySelectorAll(".mode-extra").forEach(function (box) {
    box.classList.toggle("show", checked && box.dataset.mode === checked.value);
  });
}
document.querySelectorAll("input[name='cutoff_mode']").forEach(function (radio) {
  radio.addEventListener("change", updateCutoffBoxes);
});
updateCutoffBoxes();


// ---- 5. Spinner on the Analyze button while the server works ----
var form = document.getElementById("analyze");
if (form) {
  form.addEventListener("submit", function () {
    var button = document.getElementById("submitBtn");
    button.classList.add("loading");
    button.querySelector(".btn-label").textContent = "Analyzing…";
  });
}


// ---- 6. Pick which sample to show in the gene-by-gene breakdown ----
var picker = document.getElementById("samplePicker");
if (picker) {
  picker.addEventListener("change", function () {
    document.querySelectorAll(".breakdown").forEach(function (box) {
      box.classList.toggle("active", box.dataset.index === picker.value);
    });
  });
}


// ---- 7. Download the results table as a CSV file ----
var downloadButton = document.getElementById("downloadCsv");
if (downloadButton) {
  downloadButton.addEventListener("click", function () {
    var lines = [];
    document.querySelectorAll("#resultsTable tr").forEach(function (row) {
      var cells = [];
      row.querySelectorAll("th, td").forEach(function (cell) {
        var text = cell.textContent.trim().replace(/"/g, '""');
        cells.push('"' + text + '"');
      });
      lines.push(cells.join(","));
    });
    var blob = new Blob([lines.join("\n")], { type: "text/csv" });
    var link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "panel_results.csv";
    link.click();
  });
}
