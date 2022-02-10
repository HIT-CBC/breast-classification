#!/usr/bin/env python
# coding=utf-8
"""
Transforms for images
"""
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, auc, precision_recall_curve
from sklearn.preprocessing import label_binarize
import torch
from torchvision.transforms import functional as F
from torchvision.utils import make_grid
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import copy
import math
import sys
import random
import os

class EarlyStopping:

    def __init__(self, patience=20, stop_epoch=80, verbose=False):
        self.patience = patience
        self.stop_epoch = stop_epoch
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf

    def __call__(self, epoch, val_loss, model, ckpt_name="checkpoint.pt"):
        score = -val_loss
        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model, ckpt_name)
        elif score < self.best_score:
            self.counter += 1
            print(f"EarlyStopping counter: {self.counter} out of {self.patience}")
            if self.counter >= self.patience and epoch > self.stop_epoch:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model, ckpt_name)
            self.counter = 0

    def save_checkpoint(self, val_loss, model, ckpt_name):
        if self.verbose:
            print(
                f"Validation loss decreased ({self.val_loss_min:.6f}-->{val_loss:.6f}). Saving model ..."
            )
        torch.save(model.state_dict(), ckpt_name)
        self.val_loss_min = val_loss

def imshow(inp, save_name):
    inp = inp.numpy().transpose((1, 2, 0))
    mean = np.array([0.5, 0.5, 0.5])
    std = np.array([0.1, 0.1, 0.1])
    inp = std * inp + mean
    inp = np.clip(inp, 0, 1)
    plt.imshow(inp)
    plt.savefig(save_name)


class ToTensor(object):
    """
    pil list to tensor list
    """
    def __call__(self, pil_list):
        images = [F.to_tensor(x) for x in pil_list]
        return torch.stack(images)


class Normalize(object):
    """Normalize"""
    def __call__(self, pil_list):
        images = [F.normalize(x, (0.5, 0.5, 0.5), (0.1, 0.1, 0.1))
                  for x in pil_list]
        return torch.stack(images)


class RandomHorizontalFlip(object):
    """Horizontal flip"""
    def __init__(self, p=0.5):
        self.p = p

    def __call__(self, pil_list):
        images = []
        for x in pil_list:
            if torch.rand(1) < self.p:
                images.append(F.hflip(x))
            else:
                images.append(x)
        return torch.stack(images)


class MyRotationTrans:
    def __init__(self, angles):
        self.angles = angles

    def __call__(self, x):
        angle = random.choice(self.angles)
        return F.rotate(x, angle)


def normalize_np(np_array, normalizer):
    source_array = staintools.LuminosityStandardizer.standardize(np_array)
    transformed = normalizer.transform(source_array)
    return transformed


class Compose(object):
    """
    self defined Compose like transforms.Compose
    """

    def __init__(self, transforms):
        self.transforms = transforms

    def __call__(self, pil_list):
        for t in self.transforms:
            pil_list = t(pil_list)
        return pil_list


def warmup_lr_scheduler(optimizer, warmup_iters, warmup_factor):
    """learning rate warmup"""

    def f(x):
        if x >= warmup_iters:
            return 1
        alpha = float(x) / warmup_iters
        return warmup_factor * (1 - alpha) + alpha
    return torch.optim.lr_scheduler.LambdaLR(optimizer, f)


@torch.no_grad()
def eval_process(epoch, model, criterion, dataloader, device,early_stopping=None,results_dir=None):              
    print("Epoch %d Validation......" % epoch)
    cpu_device = torch.device("cpu")
    model.eval()
    preds = []
    labels = []
    running_loss = 0.
    for images, targets, _ in tqdm(dataloader):
        images = images.to(device)
        labels.extend(targets.numpy().tolist())
        targets = targets.to(device)
        logits=model(images)
        preds.extend(logits.cpu().numpy().tolist())
        loss = criterion(logits, targets)
        running_loss += loss.item() * images.size(0)
    this_epoch_loss = running_loss / len(dataloader.dataset)
    print("Val loss: %.4f" % this_epoch_loss)
    n_classes=3
    preds = np.vstack(preds)
    labels = label_binarize(labels, classes=[0, 1, 2])
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(labels[:, i], preds[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    print(roc_auc)
    if early_stopping:
        assert results_dir
        early_stopping(epoch, this_epoch_loss, model,
                       ckpt_name=os.path.join(results_dir, "checkpoint.pt"))

        if early_stopping.early_stop:
            print("Early stopping")
            return True, roc_auc,this_epoch_loss

    return False, roc_auc,this_epoch_loss




@torch.no_grad()
def test_process(model, criterion, dataloader, device):
    cpu_device = torch.device("cpu")
    model.eval()
    preds = []
    preds_probs = []
    labels = []
    all_case_names = []
    running_loss = 0.
    for images, targets, case_names in tqdm(dataloader):
      all_case_names.extend(case_names)
      images = images.to(device)
      labels.extend(targets.numpy().tolist())
      targets = targets.to(device)
      logits = model(images)
      pred_probs = torch.softmax(logits, 1)
      preds_probs.extend(pred_probs.cpu().numpy().tolist())
      preds.extend(logits.cpu().numpy().tolist())
      loss = criterion(logits, targets)
      running_loss += loss.item() * images.size(0)
      this_epoch_loss = running_loss / len(dataloader.dataset)
      print("test loss: %.4f" % this_epoch_loss)  
    n_classes=3
    preds = np.vstack(preds)
    labels_ = label_binarize(labels, classes=[0, 1, 2])
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(labels_[:, i], preds[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    print(roc_auc)
 
    df_case_names = pd.DataFrame({"case": all_case_names})
    np_preds_probs = np.vstack(preds_probs)
    df_preds_probs = pd.DataFrame(np_preds_probs)
    df_labels = pd.DataFrame({"label": labels})
    df_total = pd.concat((df_case_names, df_preds_probs, df_labels), 1)
    df_total.rename(columns={0:"prob1", 1:"prob2", 2:"prob3"}, inplace=True)
    return roc_auc,this_epoch_loss,df_total



def train_process(model, criterion, optimizer, lr_sche, dataloaders,
                  num_epochs, use_tensorboard, device,
                  save_model_path, record_iter,
                  writer=None):
    model.train()

    best_score = None
    best_state_dict = copy.deepcopy(model.state_dict())
    early_stopping = EarlyStopping(patience=20, stop_epoch=80, verbose=True)

    for epoch in range(num_epochs):
        lr_scheduler = None
        running_loss = 0.0
        print("====Epoch{0}====".format(epoch))
        if epoch == 0:
            warmup_factor = 1. / 1000
            warmup_iters = min(1000, len(dataloaders["train"]) - 1)
            lr_scheduler = warmup_lr_scheduler(
                optimizer, warmup_iters, warmup_factor
            )

        for i, (images, targets, _) in enumerate(tqdm(dataloaders["train"])):
            images = images.to(device)
            targets = targets.to(device)
            optimizer.zero_grad()
         
            logits=model(images)
            loss = criterion(logits, targets)

            if not math.isfinite(loss.item()):
                print("Loss is {}, stopping training".format(loss.item()))
                sys.exit(1)

            loss.backward()
            optimizer.step()
            if lr_scheduler is not None:
                lr_scheduler.step()

            running_loss += loss.item() * images.size(0)

            lr = optimizer.param_groups[0]["lr"]

            if (i + 1) % record_iter == 0:
                to_date_cases = (i + 1) * images.size(0)
                tmp_loss = running_loss / to_date_cases
                print("Epoch{0} loss:{1:.4f}".format(epoch, tmp_loss))
                
                if use_tensorboard:
                    writer.add_scalar("Train loss",
                                      tmp_loss,
                                      epoch * len(dataloaders["train"]) + i)
                    writer.add_scalar("lr", lr,
                                      epoch * len(dataloaders["train"]) + i)

        stop, val_auc, val_loss = eval_process(
            epoch, model, criterion, dataloaders["val"],
            device, early_stopping,save_model_path
        )

        if use_tensorboard:
            writer.add_scalar(
                "validataion AUC", val_auc, global_step=epoch
            )
            writer.add_scalar(
                "validation loss", val_loss, global_step=epoch
            )

        model.train()
        if stop:
            break
    print("Training Done!")
    

def calculate_objective(pred, target):
    target = target.float()
    pred = torch.clamp(pred, min=1e-5, max=1. - 1e-5).squeeze(-1)
    neg_log_likelihood = -1. * (target * torch.log(pred) + (1. - target) * torch.log(1. - pred))
    neg_log_likelihood = neg_log_likelihood.mean()
    return neg_log_likelihood

