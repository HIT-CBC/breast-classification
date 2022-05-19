# A multimodal model for the diagnosis and subgroups in luminal breast cancer  

We developed a sequential multimodal model for the diagnosis and subgroups of luminal BC. A weakly supervised deep learning framework named BSNet was trained for     diagnosis in luminal BC based on multi-view unlabeled mammography of 2321 mammography cases (9284 images) and validated on the real-world external heterogeneous cohor


## Getting Started
Python3, pytorch>=1.8.0,torchvision>=0.7.0 are required for the current codebase
```
pip3 install torch==1.8.2+cu102 torchvision==0.9.2+cu102 torchaudio===0.8.2 -f https://download.pytorch.org/whl/lts/1.8/torch_lts.html
```
## Data preparation

The mammography images are placed according to the following categories
```
------Patient 1
         ------Patient 1_L_CC.png
         ------Patient 1_R_CC.png
         ------Patient 1_L_MLO.png
         ------Patient 1_R_MLO.png
------Patient 2
         ------Patient 2_L_CC.png
         ------Patient 2_R_CC.png
         ------Patient 2_L_MLO.png
         ------Patient 2_R_MLO.png
         
------Patient 3
         ------Patient 3_L_CC.png
         ------Patient 3_R_CC.png
         ------Patient 3_L_MLO.png
         ------Patient 3_R_MLO.png
...         

```

## Demo


* Train model:`bash ./run.sh`
* validate BSNet model:  `python test.py`  

## R code

Folders 1, 2, 3, and 4 are R code and data

* Folder 1 is the data preprocessing code and related data

* Folder 2 is the feature selection code and related data

* Folder 3 is the code and related data for the module selection

* Folder 4 is the clustering and typing code and related data
