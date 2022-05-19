#!/bin/bash


LOG="/home/mengyan/data/ori/resnet22dataset/fold_0/log"
NUM_EPOCHS=200
SAVE_MODEL_PATH="/home/mengyan/data/ori/resnet22dataset/fold_0"
NUM_WORKERS=10
BATCH_SIZE=16

TRAIN_CSV='/home/mengyan/data/ori/resnet22dataset/fold_0/Train.csv'
VAL_CSV='/home/mengyan/data/ori/resnet22dataset/fold_0/Val.csv'
TEST_CSV='/home/mengyan/data/ori/test.csv'

python run.py \
    --train_csv_path $TRAIN_CSV \
    --val_csv_path $VAL_CSV \
    --test_csv_path $TEST_CSV \
    --logdir $LOG \
    --num_epochs $NUM_EPOCHS \
    --save_model_path $SAVE_MODEL_PATH \
    --num_workers $NUM_WORKERS \
    --batch_size $BATCH_SIZE \


