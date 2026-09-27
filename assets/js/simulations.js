document.querySelectorAll('.simulation-viewer').forEach(function (viewer) {
  viewer.addEventListener('toggle', function () {
    var frame = viewer.querySelector('iframe');
    if (viewer.open && !frame.getAttribute('src')) {
      frame.setAttribute('src', frame.dataset.src);
    }
  });
});

document.querySelectorAll('.simulations-page video').forEach(function (video) {
  video.addEventListener('play', function () {
    document.querySelectorAll('.simulations-page video').forEach(function (other) {
      if (other !== video) other.pause();
    });
  });
});
