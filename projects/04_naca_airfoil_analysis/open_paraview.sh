#!/bin/bash
# =========================================================================
#  NACA 0012 — unified ParaView launcher
#
#  Usage:
#    ./open_paraview.sh                        Steady case, interactive (default)
#    ./open_paraview.sh steady                 Steady case, interactive
#    ./open_paraview.sh steady --screenshot s  Steady case, save screenshot
#    ./open_paraview.sh transient              Transient case, interactive
#    ./open_paraview.sh transient render       Batch render frames + MP4
#    ./open_paraview.sh transient render-only  Frames only
#    ./open_paraview.sh transient video-only   ffmpeg only (existing frames)
# =========================================================================
set -euo pipefail

cd "${0%/*}" || exit 1

CASE="${1:-steady}"
shift 2>/dev/null || true

# --- Transient animation tunables ----------------------------------------
FRAME_DIR="transient_animation/frames"
VIDEO="transient_animation/naca_vorticity_animation.mp4"
FPS=24
RESOLUTION="1920 1080"
SKIP_BEFORE="2.0"
VORT_RANGE="-30 30"
ZOOM="1.5"
CAMERA_X="2.0"
CAMERA_Y="-0.1"

usage() {
    cat <<EOF
Usage: $(basename "$0") [CASE] [MODE] [extra args...]

Cases:
  steady      Steady-state RANS visualisation (default)
  transient   Transient vorticity animation

Transient modes (second argument):
  (empty)       Open ParaView GUI interactively (default)
  render        Export PNGs with pvbatch, then assemble MP4
  render-only   Export PNGs only
  video-only    Assemble MP4 from existing frames

Extra arguments are forwarded to the Python script.
Tunable knobs (resolution, fps, vort range, etc.) are at the top of this script.
EOF
}

check_cmd() {
    if ! command -v "$1" &>/dev/null; then
        echo "ERROR: '$1' not found in PATH." >&2
        return 1
    fi
}

# --- Steady ---------------------------------------------------------------

run_steady() {
    touch airfoil.foam
    check_cmd paraview
    echo "Opening steady-state case in ParaView..."
    paraview --script paraview_airfoil.py
}

# --- Transient helpers ----------------------------------------------------

render_frames() {
    check_cmd pvbatch
    echo "========================================"
    echo "  Exporting animation frames (pvbatch)"
    echo "========================================"
    # shellcheck disable=SC2086
    pvbatch transient_animation/paraview_streamlines.py \
        --case transient_animation/airfoil.foam \
        --batch \
        --output-dir "$FRAME_DIR" \
        --resolution $RESOLUTION \
        --skip-before "$SKIP_BEFORE" \
        --vort-range $VORT_RANGE \
        --zoom "$ZOOM" \
        --camera-x "$CAMERA_X" \
        --camera-y "$CAMERA_Y" \
        --streamlines \
        "$@"
    echo "Frames saved to $FRAME_DIR/"
}

make_video() {
    check_cmd ffmpeg
    local nframes
    nframes=$(find "$FRAME_DIR" -maxdepth 1 -name 'frame_*.png' 2>/dev/null | wc -l)
    if [ "$nframes" -eq 0 ]; then
        echo "ERROR: no frame_*.png files in $FRAME_DIR/" >&2
        exit 1
    fi
    echo "========================================"
    echo "  Assembling $nframes frames → $VIDEO (${FPS} fps)"
    echo "========================================"
    ffmpeg -y -framerate "$FPS" \
        -i "$FRAME_DIR/frame_%05d.png" \
        -c:v libx264 -pix_fmt yuv420p \
        -crf 18 -preset slow \
        -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
        "$VIDEO"
    echo "Video saved: $VIDEO"
}

run_transient() {
    touch transient_animation/airfoil.foam
    local MODE="${1:-interactive}"
    shift 2>/dev/null || true

    case "$MODE" in
        interactive)
            check_cmd paraview
            echo "Opening transient case in ParaView..."
            paraview --script transient_animation/paraview_streamlines.py
            ;;
        render)
            render_frames "$@"
            make_video
            ;;
        render-only)
            render_frames "$@"
            ;;
        video-only)
            make_video
            ;;
        *)
            echo "Unknown transient mode: $MODE" >&2
            usage
            exit 1
            ;;
    esac
}

# --- Dispatch -------------------------------------------------------------

case "$CASE" in
    steady)                 run_steady "$@" ;;
    transient)              run_transient "$@" ;;
    -h|--help|help)         usage ;;
    *)
        echo "Unknown case: $CASE" >&2
        usage
        exit 1
        ;;
esac
