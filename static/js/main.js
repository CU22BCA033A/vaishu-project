// main.js - a few small extras. Every page also works without JavaScript.


// 1. Light / Dark switch in the header (the choice is remembered in the browser)
function currentTheme() {
  var chosen = document.documentElement.getAttribute("data-theme");
  if (chosen) return chosen;
  return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
}
function markThemeButtons() {
  var buttons = document.querySelectorAll("[data-theme-choice]");
  for (var i = 0; i < buttons.length; i++) {
    var on = buttons[i].getAttribute("data-theme-choice") === currentTheme();
    buttons[i].setAttribute("aria-pressed", on ? "true" : "false");
  }
}
var themeButtons = document.querySelectorAll("[data-theme-choice]");
for (var t = 0; t < themeButtons.length; t++) {
  themeButtons[t].addEventListener("click", function () {
    var theme = this.getAttribute("data-theme-choice");
    document.documentElement.setAttribute("data-theme", theme);
    try { localStorage.setItem("theme", theme); } catch (e) {}
    markThemeButtons();
  });
}
markThemeButtons();


// 1b. The drop zone: show the chosen file name, highlight while dragging
var drop = document.getElementById("drop");
if (drop) {
  var fileInput = drop.querySelector("input[type='file']");
  var nameLabel = drop.querySelector(".file-name");
  fileInput.addEventListener("change", function () {
    if (fileInput.files.length > 0) {
      var file = fileInput.files[0];
      nameLabel.textContent = file.name + " · " + Math.max(1, Math.round(file.size / 1024)) + " KB";
      drop.classList.add("has-file");
    }
  });
  drop.addEventListener("dragover", function () { drop.classList.add("over"); });
  drop.addEventListener("dragleave", function () { drop.classList.remove("over"); });
  drop.addEventListener("drop", function () { drop.classList.remove("over"); });
}


// 2. Only show the extra boxes for the cutoff option that is selected
function showCutoffExtras() {
  var picked = document.querySelector("input[name='cutoff_mode']:checked");
  var boxes = document.querySelectorAll(".extra");
  for (var i = 0; i < boxes.length; i++) {
    var show = picked && boxes[i].getAttribute("data-for") === picked.value;
    boxes[i].classList.toggle("show", show);
  }
}
var cutoffRadios = document.querySelectorAll("input[name='cutoff_mode']");
for (var i = 0; i < cutoffRadios.length; i++) {
  cutoffRadios[i].addEventListener("change", showCutoffExtras);
}
showCutoffExtras();


// 3. Disable the button while the server is working, so it is not clicked twice
var form = document.getElementById("analyzeForm");
if (form) {
  form.addEventListener("submit", function () {
    var button = document.getElementById("submitBtn");
    button.disabled = true;
    button.textContent = "Analyzing…";
  });
}


// 4. Choose which sample the gene-by-gene bars show
var picker = document.getElementById("samplePicker");
if (picker) {
  picker.addEventListener("change", function () {
    var boxes = document.querySelectorAll(".breakdown");
    for (var i = 0; i < boxes.length; i++) {
      boxes[i].classList.toggle("show", boxes[i].getAttribute("data-index") === picker.value);
    }
  });
}


// 5. Download the results table as a CSV file
var downloadBtn = document.getElementById("downloadBtn");
if (downloadBtn) {
  downloadBtn.addEventListener("click", function () {
    var lines = [];
    var rows = document.querySelectorAll("#resultsTable tr");
    for (var r = 0; r < rows.length; r++) {
      var cells = rows[r].querySelectorAll("th, td");
      var values = [];
      for (var c = 0; c < cells.length; c++) {
        values.push('"' + cells[c].textContent.trim().replace(/"/g, '""') + '"');
      }
      lines.push(values.join(","));
    }
    var file = new Blob([lines.join("\n")], { type: "text/csv" });
    var link = document.createElement("a");
    link.href = URL.createObjectURL(file);
    link.download = "panel_results.csv";
    link.click();
  });
}
