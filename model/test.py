#!/usr/bin/env python
# coding=utf-8
from torch.utils.data import DataLoader
from torchvision import transforms
from newdata import MILBreastDataset
from utils import test_process
from model_res22 import ResNet22
import torch
import torch.nn as nn


device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")

transform = transforms.Compose([
    transforms.Resize((512, 512)),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.5], std=[0.5])
])
print("====Preparing Model====")
model = ResNet22(batch_size=50, dropout=True)
model = model.to(device)
if torch.cuda.device_count() > 1:
    model = nn.DataParallel(model)
model.load_state_dict(torch.load("/home/mengyan/data/ori/resnet22dataset/fold_0/checkpoint.pt"))

print("====Dataset====")
test_dataset = MILBreastDataset("/home/mengyan/data/exvali.csv",
                                 transform=transform)
test_loader = DataLoader(test_dataset, batch_size=50,
                         shuffle=False, num_workers=12)

criterion = torch.nn.CrossEntropyLoss()
print("====Testing====")
test_auc, test_loss, df = test_process(
    model, criterion, test_loader, device
)
df.to_csv("/home/mengyan/data/valiresnet22_4_prediction.csv", index=False)
print("Done")

