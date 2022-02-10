#!/usr/bin/env python
# coding=utf-8
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, WeightedRandomSampler
from torchvision import transforms
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm

import argparse
import os
import warnings
import pickle

warnings.filterwarnings("ignore")

from newdata import MILBreastDataset, make_weights_for_balanced_classes_split
from utils import train_process, MyRotationTrans, grid_show #Cutout, PixelReg
#from model import AttBMNet,MAX, ResNet22
from model_res22 import MAX
from utils import Compose, ToTensor

device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
parser = argparse.ArgumentParser(description="Running script")

parser.add_argument(
    "--train_csv_path", type=str, default="./train.csv",
    help="train case path csv, default ./train.csv"
)
parser.add_argument(
    "--val_csv_path", type=str, default="./val.csv",
    help="val case path csv, default ./val.csv"
)
parser.add_argument(
    "--test_csv_path", type=str, default="./test.csv",
    help="test case path csv, default ./test.csv"
)
#parser.add_argument(
    #"--lmdb_dir", type=str, default="./lmdb_dir",
    #help="lmdb directory, default ./lmdb_dir"
#)
parser.add_argument(
    "--batch_size", type=int, default=4,
    help="batch size, default 16"
)
parser.add_argument(
    "--num_workers", type=int, default=4,
    help="num workers, default 16"
)
#parser.add_argument(
    #"--k", type=int, default=40,
    #help="random sample k patches from each slide, default 40"
#)
parser.add_argument(
    "--lr", type=float, default=0.0005,
    help="learning rate, default 0.001"
)
parser.add_argument(
    "--momentum", type=float, default=0.9,
    help="SGD momentum, default 0.9"
)
parser.add_argument(
    "--weight_decay", type=float, default=5e-4,
    help="SGD weight decay, default 5e-4"
)
parser.add_argument(
    "--gamma", type=float, default=0.1,
    help="StepLR gamma value, default 0.1"
)
parser.add_argument(
    "--num_epochs", type=int, default=10,
    help="num of epochs, default 500"
)
parser.add_argument(
    "--save_model_path", type=str,
    default="/home/dl/code/MengYan/dataloader_exercise/best_model.pth",
    help="model saving path"
)
parser.add_argument(
    "--record_iter", type=int, default=10,
    help="print record iter, default 10"
)

parser.add_argument(
    "--logdir", type=str, default="./logs",
    help="tensorboard log dir, default ./logs"
)

def main(args):
    
    rotation = MyRotationTrans([0, 90, 180, 270])
    transform = transforms.Compose([
        transforms.Resize((512, 512)),
        transforms.RandomHorizontalFlip(p=0.5),
        rotation,
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.5], std=[0.5]),
    ])
    val_transform = transforms.Compose([
        transforms.Resize((512, 512)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.5], std=[0.5])
    ])
    print("========Preparing Dataset========")
    train_dataset = MILBreastDataset(args.train_csv_path, transform)
    weights = make_weights_for_balanced_classes_split(train_dataset)
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size,
                              num_workers=args.num_workers, sampler=WeightedRandomSampler(weights, len(weights)))
    val_dataset = MILBreastDataset(args.val_csv_path, transform=val_transform)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size,
                            shuffle=False, num_workers=args.num_workers)
    test_dataset = MILBreastDataset(args.test_csv_path,  transform=val_transform)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size,
                             shuffle=False, num_workers=args.num_workers)
    dataloaders = {"train": train_loader, "val": val_loader, "test": test_loader}
    print("========Dataset Done========")

    print("========Preparing Model========")
    model = MAX(batch_size=args.batch_size,dropout=True)
    model = model.to(device)
    if torch.cuda.device_count() > 1:
        model = nn.DataParallel(model)
    print("========Model Done========")

    params = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.SGD(params, lr=args.lr,
                                momentum=args.momentum,
                                weight_decay=args.weight_decay)
    
    lr_scheduler = None
    optimizer = torch.optim.Adam(params, lr=1e-4, weight_decay=5e-4)
    loss_weights = torch.tensor([0.2, 0.25, 0.55])
    loss_weights = loss_weights.to(device)
    criterion = nn.CrossEntropyLoss(weight=loss_weights)
    
    print("========Start Training========")
    logdir = args.logdir
    os.makedirs(logdir, exist_ok=True)
    writer = SummaryWriter(logdir)
    train_process(model=model, criterion=criterion, optimizer=optimizer,
                  lr_sche=lr_scheduler, dataloaders=dataloaders, writer=writer,
                  num_epochs=args.num_epochs, use_tensorboard=args.use_tensorboard,
                  device=device, save_model_path=args.save_model_path,
                  record_iter=args.record_iter)


if __name__ == "__main__":
    main(args)

