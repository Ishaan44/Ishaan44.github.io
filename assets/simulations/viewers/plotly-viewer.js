(() => {
  "use strict";
  const figure = JSON.parse(document.getElementById("figure-data").textContent);
  const plot = document.getElementById("plot");
  const status = document.getElementById("plot-status");
  const layout = { ...figure.layout };
  // Keep the archived figure intact; move only its presentation outside the plot.
  delete layout.width;
  delete layout.height;
  delete layout.title;
  delete layout.annotations;
  layout.autosize = true;
  layout.showlegend = false;
  layout.margin = { l: 12, r: 12, t: 30, b: 32, pad: 0 };
  const fitScale = () => Math.max(1.08, plot.clientHeight / Math.max(1, plot.clientWidth) * 1.1);
  let cameraScale = fitScale();
  const initialCamera = figure.layout.scene.camera;
  layout.scene = {
    ...figure.layout.scene,
    camera: {
      ...initialCamera,
      eye: Object.fromEntries(Object.entries(initialCamera.eye).map(([axis, value]) => [axis, value * cameraScale])),
    },
  };
  const config = { ...figure.config, responsive: true, displaylogo: false };

  Plotly.newPlot(plot, figure.data, layout, config).then(() => {
    status.textContent = "";
    const controls = document.getElementById("trace-controls");
    const options = document.getElementById("trace-options");
    figure.data.forEach((trace, index) => {
      if (trace.showlegend === false || !trace.name) return;
      const label = document.createElement("label");
      label.className = "trace-option";
      const input = document.createElement("input");
      input.type = "checkbox";
      input.checked = trace.visible !== "legendonly" && trace.visible !== false;
      const swatch = document.createElement("span");
      swatch.className = "trace-swatch";
      swatch.setAttribute("aria-hidden", "true");
      const color = trace.line?.color || trace.marker?.color;
      swatch.style.backgroundColor = typeof color === "string" ? color : "#555";
      const text = document.createElement("span");
      text.className = "trace-label";
      text.textContent = trace.name;
      input.addEventListener("change", () => {
        const indices = trace.legendgroup
          ? figure.data.map((item, itemIndex) => item.legendgroup === trace.legendgroup ? itemIndex : -1).filter(itemIndex => itemIndex >= 0)
          : [index];
        Plotly.restyle(plot, { visible: input.checked ? true : "legendonly" }, indices);
      });
      label.append(input, swatch, text);
      options.append(label);
    });
    controls.hidden = !options.childElementCount;
    const observer = new ResizeObserver(() => {
      const nextScale = fitScale();
      if (Math.abs(nextScale - cameraScale) > 0.005) {
        const eye = plot.layout.scene.camera.eye;
        const adjusted = Object.fromEntries(Object.entries(eye).map(([axis, value]) => [axis, value * nextScale / cameraScale]));
        cameraScale = nextScale;
        Plotly.relayout(plot, { "scene.camera.eye": adjusted });
      }
      Plotly.Plots.resize(plot);
    });
    observer.observe(plot);
  }).catch(error => {
    status.textContent = "The saved figure could not be displayed. Reload this page to try again.";
    console.error(error);
  });
})();
