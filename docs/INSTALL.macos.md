# HOWARD Installation guide on MacOS

Follow these steps to install and set up HOWARD on macOS.

## Step 1: Create a Conda Environment

First, create a dedicated Conda environment for HOWARD:

```bash
conda create --name howard python=3.10
conda activate howard
```

This ensures that HOWARD and its dependencies are isolated from other projects.

## Step 2: Install OpenMP and LLVM

Install the required OpenMP and LLVM libraries using Homebrew:

```bash
brew install llvm
brew install libomp
```

These libraries are necessary for compiling and running HOWARD efficiently.

## Step 3: Define Environment Variables

Set the required environment variables to link the installed libraries:

```bash
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
```

These flags ensure that the compiler can locate the OpenMP library and its headers.

## Step 4: Install HOWARD Using pip

Install HOWARD in editable mode using pip:

```bash
python -m pip install -e .
```

## Step 5: Verify the Installation

Test the HOWARD installation to ensure everything is set up correctly:

```bash
howard --help
```

If the installation is successful, this command will display the HOWARD help message.
