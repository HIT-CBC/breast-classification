#!/usr/bin/env python
# coding=utf-8

import torch
from torch.utils.data import Dataset
from torchvision import transforms
import pandas as pd
from PIL import Image
from tqdm import tqdm
import os
from glob import glob
from torch.utils.data import DataLoader

class MILBreastDataset(Dataset):

    def __init__(self, csv_path, transform=None):
      
        df = pd.read_csv(csv_path,encoding="utf_8_sig")
        self.case_dirs = list(df.case)
        self.labels = list(df.label)
        self.transform = transform
        self.slide_cls_ids = {0: self.labels.count(0),
                              1: self.labels.count(1),
                              2: self.labels.count(2)}

    def __len__(self):
      
        return len(self.case_dirs)

    def __getitem__(self, index):
        case_dir = self.case_dirs[index]
        img_paths_of_this_case = glob(case_dir + "/*.png")
        samples=[]
        for img_path in img_paths_of_this_case:
            img = Image.open(img_path)
            samples.append(img)
        if self.transform is not None:
            trans_samples = []
            for sample in samples:
                trans_samples.append(self.transform(sample))
            samples = torch.stack(trans_samples)
        return samples, self.labels[index], case_dir

    def getlabel(self, idx):
        return self.labels[idx]

def make_weights_for_balanced_classes_split(dataset):
    N = float(len(dataset))                                           
    weight_per_class = [N/(dataset.slide_cls_ids[c]) for c in range(len(dataset.slide_cls_ids))]                                                                                                     
    weight = [0] * int(N)                                           
    for idx in range(len(dataset)):   
        y = dataset.getlabel(idx)                        
        weight[idx] = weight_per_class[y]                                  
    return torch.DoubleTensor(weight)


if __name__ == "__main__":
  
    transform = transforms.Compose([
        transforms.Resize((512, 512)),
        transforms.ToTensor()
    ])
    dataset = MILBreastDataset(csv_path, transform)
    trainloader = DataLoader(dataset, batch_size=16, shuffle=True, num_workers=0)
   
