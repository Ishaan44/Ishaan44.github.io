document.querySelectorAll('.simulation-viewer').forEach(function (viewer) {
  viewer.addEventListener('toggle', function () {
    var frame = viewer.querySelector('iframe');
    if (viewer.open && !frame.getAttribute('src')) {
      frame.setAttribute('src', frame.dataset.src);
    }
  });
});

// Static PNG frames avoid the platform video decoder entirely.
document.querySelectorAll('.simulation-player').forEach(function (player) {
  var image = player.querySelector('.simulation-frame');
  var play = player.querySelector('.simulation-play');
  var icon = play.querySelector('img');
  var seek = player.querySelector('.simulation-seek');
  var time = player.querySelector('.simulation-time');
  var status = player.querySelector('.simulation-status');
  var fullscreen = player.querySelector('.simulation-fullscreen');
  var manifestUrl = new URL(player.dataset.sequence, window.location.href);
  var manifestPromise, manifest, starts, frameIndex = 0;
  var cache = new Map();
  var playing = false, generation = 0, animation = 0;
  var position = 0, anchor = null, duration = Number(player.dataset.duration) * 1000;

  function format(ms) {
    var seconds = Math.floor(ms / 1000);
    return Math.floor(seconds / 60) + ':' + String(seconds % 60).padStart(2, '0');
  }
  function displayTime(value) {
    seek.value = value / 1000;
    time.textContent = format(value) + ' / ' + format(duration);
    seek.setAttribute('aria-valuetext', format(value) + ' of ' + format(duration));
  }
  function setPlaying(value) {
    playing = value;
    player.dataset.playing = String(value);
    play.setAttribute('aria-label', value ? 'Pause' : 'Play');
    play.title = value ? 'Pause' : 'Play';
    icon.src = value ? icon.dataset.pause : icon.dataset.play;
  }
  function currentPosition() {
    if (anchor === null) return position;
    return Math.min(starts[frameIndex + 1], position + performance.now() - anchor);
  }
  function cancel() {
    generation++;
    cancelAnimationFrame(animation);
    anchor = null;
  }
  function pause() {
    position = currentPosition();
    // Finish an in-flight seek even when playback is paused while it loads.
    cancelAnimationFrame(animation);
    anchor = null;
    setPlaying(false);
    status.textContent = '';
    displayTime(position);
  }
  function loadManifest() {
    if (!manifestPromise) {
      manifestPromise = fetch(manifestUrl).then(function (response) {
        if (!response.ok) throw new Error('Frame index unavailable');
        return response.json();
      }).then(function (data) {
        if (!Array.isArray(data.frames) || !data.frames.length) throw new Error('Empty animation');
        manifest = data;
        starts = [0];
        data.frames.forEach(function (frame) { starts.push(starts[starts.length - 1] + frame.duration_ms); });
        duration = starts[starts.length - 1];
        seek.max = duration / 1000;
        return data;
      }).catch(function (error) { manifestPromise = null; throw error; });
    }
    return manifestPromise;
  }
  function loadFrame(index) {
    if (!cache.has(index)) {
      var frame = new Image();
      frame.src = new URL(manifest.frames[index].file, manifestUrl).href;
      cache.set(index, frame.decode().then(function () { return frame; }).catch(function (error) {
        cache.delete(index);
        throw error;
      }));
    }
    return cache.get(index);
  }
  function findFrame(value) {
    var low = 0, high = manifest.frames.length - 1;
    while (low < high) {
      var mid = Math.floor((low + high + 1) / 2);
      if (starts[mid] <= value) low = mid;
      else high = mid - 1;
    }
    return low;
  }
  async function showFrame(index, token) {
    var frame = await loadFrame(index);
    if (token !== generation) return false;
    // Keep the previous picture visible until the replacement has decoded.
    image.src = frame.src;
    image.dataset.frame = index;
    frameIndex = index;
    for (var next = index + 1; next < Math.min(index + 4, manifest.frames.length); next++) {
      loadFrame(next).catch(function () {});
    }
    cache.forEach(function (_, key) { if (key < index - 1 || key > index + 4) cache.delete(key); });
    return true;
  }
  function fail(token) {
    if (token !== generation) return;
    pause();
    status.textContent = 'Frame unavailable.';
  }
  function tick(token) {
    if (!playing || token !== generation) return;
    var value = currentPosition();
    displayTime(value);
    if (value < starts[frameIndex + 1]) {
      animation = requestAnimationFrame(function () { tick(token); });
      return;
    }
    position = starts[frameIndex + 1];
    anchor = null;
    if (frameIndex === manifest.frames.length - 1) { pause(); return; }
    status.textContent = 'Loading...';
    showFrame(frameIndex + 1, token).then(function (shown) {
      if (!shown) return;
      status.textContent = '';
      if (playing) {
        anchor = performance.now();
        tick(token);
      }
    }).catch(function () { fail(token); });
  }
  async function go(value, resume) {
    cancel();
    var token = generation;
    position = Math.max(0, Math.min(value, duration));
    setPlaying(resume);
    if (resume) document.dispatchEvent(new CustomEvent('simulationplay', { detail: player }));
    displayTime(position);
    status.textContent = 'Loading...';
    try {
      await loadManifest();
      if (token !== generation || !await showFrame(findFrame(position), token)) return;
      status.textContent = '';
      if (resume && playing && position < duration) {
        anchor = performance.now();
        tick(token);
      } else setPlaying(false);
    } catch (_) { fail(token); }
  }
  play.addEventListener('click', function () {
    if (playing) pause();
    else go(position >= duration ? 0 : position, true);
  });
  player.querySelector('.simulation-replay').addEventListener('click', function () { go(0, true); });
  seek.addEventListener('input', function () { go(Number(seek.value) * 1000, playing); });
  document.addEventListener('simulationplay', function (event) { if (event.detail !== player && playing) pause(); });
  document.addEventListener('visibilitychange', function () { if (document.hidden && playing) pause(); });
  player.addEventListener('simulationpause', pause);
  if (!player.requestFullscreen) fullscreen.hidden = true;
  fullscreen.addEventListener('click', function () {
    var action = document.fullscreenElement ? document.exitFullscreen() : player.requestFullscreen();
    if (action) action.catch(function () {});
  });
  setPlaying(false);
  displayTime(0);
});

var recordingSelect = document.querySelector('.simulation-recording-select');
if (recordingSelect) {
  function selectRecording() {
    document.querySelectorAll('[data-recording]').forEach(function (recording) {
      recording.querySelector('.simulation-player').dispatchEvent(new Event('simulationpause'));
      recording.hidden = recording.dataset.recording !== recordingSelect.value;
    });
    var url = new URL(window.location.href);
    url.searchParams.set('id', recordingSelect.value);
    history.replaceState(null, '', url);
  }
  var requested = new URLSearchParams(window.location.search).get('id');
  if (Array.from(recordingSelect.options).some(function (option) { return option.value === requested; })) recordingSelect.value = requested;
  recordingSelect.addEventListener('change', selectRecording);
  selectRecording();
}
