"""Procedural, royalty-free ambient background music for the reels (numpy only).

A slow Fmaj7 - Em7 - Dm7 - Cmaj7 progression: warm detuned pads, a soft bass, and
sparse electric-piano notes, all through a long convolution reverb. The seed makes
the result repeatable, so re-rendering a reel gives the same soundtrack.

    pvpython tools/ambient_music.py music.wav --seconds 60
"""
import argparse
import wave

import numpy as np

RATE = 48000
BPM = 72
# (bass, pad voicing) as MIDI notes; each chord lasts two bars.
PROGRESSION = [
    (41, (57, 60, 64, 67)),  # Fmaj7(9)
    (40, (55, 59, 62, 64)),  # Em7
    (38, (53, 57, 60, 64)),  # Dm7(9)
    (36, (52, 55, 59, 62)),  # Cmaj7(9)
]


def hz(note):
    return 440.0 * 2 ** ((note - 69) / 12)


def smooth_window(length, attack, release):
    """Raised-cosine attack and release (in samples) around a flat sustain."""
    env = np.ones(length)
    a, r = min(attack, length // 2), min(release, length // 2)
    env[:a] = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, a))
    env[length - r:] = 0.5 + 0.5 * np.cos(np.linspace(0, np.pi, r))
    return env


def add(buffer, start, signal):
    end = min(len(buffer), start + len(signal))
    if end > start:
        buffer[start:end] += signal[:end - start]


def reverb(channel, rng, seconds=3.2, wet=0.38):
    """Convolve with exponentially decaying, slightly darkened noise."""
    n = int(seconds * RATE)
    t = np.arange(n) / RATE
    ir = rng.standard_normal(n) * np.exp(-6.9 * t / seconds)  # -60 dB at `seconds`
    ir = np.convolve(ir, np.ones(12) / 12, mode="same")  # gentle low-pass
    ir[: int(0.02 * RATE)] = 0  # pre-delay
    ir /= np.sqrt(np.sum(ir ** 2))
    size = 1 << int(np.ceil(np.log2(len(channel) + n)))
    tail = np.fft.irfft(np.fft.rfft(channel, size) * np.fft.rfft(ir, size), size)[: len(channel)]
    return (1 - wet) * channel + wet * tail


def compose(seconds, seed=11):
    """Stereo float array of shape (samples, 2), peak-normalised to 0.9."""
    rng = np.random.default_rng(seed)
    total = int((seconds + 4) * RATE)  # extra tail for the reverb
    left, right = np.zeros(total), np.zeros(total)
    beat = 60.0 / BPM
    chord_len = 8 * beat
    n_chords = int(np.ceil(seconds / chord_len)) + 1

    for index in range(n_chords):
        bass, voicing = PROGRESSION[index % len(PROGRESSION)]
        start = int(index * chord_len * RATE)

        # Pads: two detuned partial stacks per note, panned apart, with slow overlap.
        length = int((chord_len + 2.5) * RATE)
        t = np.arange(length) / RATE
        env = smooth_window(length, int(2.2 * RATE), int(2.8 * RATE))
        swell = 1 + 0.15 * np.sin(2 * np.pi * 0.07 * t + index)
        for note in voicing:
            f = hz(note)
            for detune, channel in ((-0.0018, left), (0.0018, right)):
                phase = rng.uniform(0, 2 * np.pi)
                wave_ = (np.sin(2 * np.pi * f * (1 + detune) * t + phase)
                         + 0.18 * np.sin(4 * np.pi * f * (1 + detune) * t + phase)
                         + 0.05 * np.sin(6 * np.pi * f * (1 + detune) * t))
                add(channel, start, 0.035 * env * swell * wave_)

        # Bass: a soft sine on beats 1 and 5 of each chord.
        for hit, level in ((0, 1.0), (4, 0.6)):
            blen = int(3.2 * beat * RATE)
            tb = np.arange(blen) / RATE
            benv = (1 - np.exp(-tb / 0.03)) * np.exp(-tb / 2.2) * smooth_window(blen, 1, int(0.4 * RATE))
            sig = 0.11 * level * benv * (np.sin(2 * np.pi * hz(bass) * tb)
                                         + 0.12 * np.sin(4 * np.pi * hz(bass) * tb))
            add(left, start + int(hit * beat * RATE), sig)
            add(right, start + int(hit * beat * RATE), sig)

        # Electric piano: sparse chord tones on eighth notes, never in the first bar.
        for eighth in range(16):
            when = index * chord_len + eighth * beat / 2
            if when < 4 * beat or when > seconds - 2 or rng.random() > 0.3:
                continue
            note = int(rng.choice(voicing)) + int(rng.choice((12, 12, 24)))
            note = min(note, 88)
            f = hz(note)
            klen = int(3.0 * RATE)
            tk = np.arange(klen) / RATE
            kenv = ((1 - np.exp(-tk / 0.004)) * np.exp(-tk / 1.1)
                    * smooth_window(klen, 1, int(0.6 * RATE)))  # no click at the cut
            tone = (np.sin(2 * np.pi * f * tk)
                    + 0.35 * np.exp(-tk / 0.25) * np.sin(4 * np.pi * f * tk)
                    + 0.08 * np.exp(-tk / 0.1) * np.sin(6 * np.pi * f * tk))
            sig = 0.07 * rng.uniform(0.55, 1.0) * kenv * tone
            pan = rng.uniform(0.25, 0.75)
            at = int(when * RATE)
            add(left, at, sig * np.sqrt(1 - pan))
            add(right, at, sig * np.sqrt(pan))

    stereo = np.stack([reverb(left, rng), reverb(right, rng)], axis=1)
    stereo = stereo[: int(seconds * RATE)]
    return 0.9 * stereo / max(np.max(np.abs(stereo)), 1e-9)


def write_wav(path, stereo):
    data = (np.clip(stereo, -1, 1) * 32767).astype("<i2")
    with wave.open(str(path), "wb") as out:
        out.setnchannels(2)
        out.setsampwidth(2)
        out.setframerate(RATE)
        out.writeframes(data.tobytes())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("output")
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()
    write_wav(args.output, compose(args.seconds, args.seed))
