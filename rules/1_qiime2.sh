#!/usr/bin/env bash
set -euo pipefail

# ---- 설정(환경변수로 덮어쓰기 가능) ----
ENV_NAME="${ENV_NAME:-qiime2}"
YML_FILE="${YML_FILE:-qiime2-amplicon-ubuntu-latest-conda.yml}"
CONDA_BIN="${CONDA_BIN:-conda}"

echo "[env] CONDA_BIN=$CONDA_BIN  ENV_NAME=$ENV_NAME  YML_FILE=$YML_FILE"

# env 존재 여부에 따라 create/update
if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[env] '$ENV_NAME' exists → update from $YML_FILE"
  "$CONDA_BIN" env update -n "$ENV_NAME" -f "$YML_FILE"
else
  echo "[env] create '$ENV_NAME' from $YML_FILE"
  "$CONDA_BIN" env create -n "$ENV_NAME" -f "$YML_FILE"
fi

# 권장 옵션(실패해도 무시)
"$CONDA_BIN" config --set channel_priority flexible || true

echo "[env] done"

