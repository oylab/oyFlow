# oyFlow
Code for analysis of flow cytometry data. Based heavily on [FlowCytometryTools](https://eyurtsev.github.io/FlowCytometryTools/):


# Install 

```
conda create --name oyflow python=3.9 scipy

conda activate oyflow
```
## if you have cuda 11.2:

`pip install git+https://github.com/oylab/oyFlow.git#egg=oyFlow`

### **add kernel to jupyter:**

`python -m ipykernel install --user --name=oyflow`


## on Linux/Mac from source:
```
git clone https://github.com/oylab/oyFlow

cd oyFlow

conda create --name oyflow python=3.9 scipy

conda activate oyflow
```
`pip install -e .`

### **add kernel to jupyter:**

`python -m ipykernel install --user --name=oyflow`
