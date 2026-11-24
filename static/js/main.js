// main.js - a few small extras. Every page also works without JavaScript.


// 1. Dark / light theme button (the choice is remembered in the browser)
var themeBtn = document.getElementById("themeBtn");
if (themeBtn) {
  themeBtn.addEventListener("click", function () {
    var now = document.documentElement.getAttribute("data-theme");
    if (!now) {
      now = window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
    }
    var next = (now === "dark") ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    try { localStorage.setItem("theme", next); } catch (e) {}
  });
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
