#!/usr/bin/env bash


  # 5회 반복
  for i in $(seq 2 2 ); do
    echo "=== Fixed Simulation repeat ${i} ==="

    # 1) 시뮬레이션 실행
    ./ns3 run "scratch/lstm_trajectory_estimation_scenario_train.cc"

    # 2) txt 파일 압축
    ARCHIVE_NAME="data_lstm_${i}.tar.gz"
    echo "Archiving txt files ⇒ ${ARCHIVE_NAME}"
    find . -maxdepth 1 -type f -name '*.txt' \
         ! -name 'CMakeLists.txt' \
         -print0 \
      | tar --null -czvf "${ARCHIVE_NAME}" --files-from -

    # 3) 원본 txt 파일 삭제
    find . -maxdepth 1 -type f -name '*.txt' \
         ! -name 'CMakeLists.txt' \
         -delete
  done

echo "✅ All fixed-run simulations completed!"
