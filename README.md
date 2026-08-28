# BGO Testing Reconstruction and Analysis

This repository contains the code used to reconstruct and analyse the BGO crystal test data.

## Prerequisites

You need:

- ROOT installed and working
- access to the raw data
- `hadd` available from your ROOT installation

Check that ROOT works before doing anything else:

```bash
root --version
```

and check where `hadd` is installed:

```bash
which hadd
```

---

# 1. Initial setup

Before running any reconstruction, execute:

```bash
./setup.sh
```

This creates the folders required by the analysis.

You only need to do this once after cloning the repository.

If the script is not executable:

```bash
chmod +x setup.sh
./setup.sh
```

The analysis expects the following folders to exist:

```text
data/
calib/
out/
pdf/
```

---

# 2. Getting the raw data

The raw data are stored on the DAQ machine in SGM.

If you are connected to the INFN network, either directly or through the VPN, copy a run with:

```bash
scp streamdaq@193.206.147.141:/home/streamdaq/Documents/BGO_test/DAQ/RUN/RAW/SDataR_RUN.root ./data/
```

Replace `RUN` with the actual run number.

For example:

```bash
scp streamdaq@193.206.147.141:/home/streamdaq/Documents/BGO_test/DAQ/260824/RAW/SDataR_260824.root ./data/
```

If you are not connected to the INFN network or VPN, you will need someone with access to copy the data for you.

## Rename the files

It is strongly recommended to remove the CoMPASS prefix from the filename.

For example:

```text
SDataR_260824.root
```

should become:

```text
260824.root
```

You can do this with:

```bash
mv data/SDataR_260824.root data/260824.root
```

The reconstruction scripts assume this simpler naming convention.

---

# 3. Keep `run_list.dat` updated

The file:

```text
run_list.dat
```

contains the association between crystal batches and DAQ runs.

Example:

```text
3
260824
260825
260826
2
260814
260815
```

This means that:

```text
Batch 3:
    260824
    260825
    260826

Batch 2:
    260814
    260815
```

Whenever a new run is taken, make sure that `run_list.dat` is updated correctly.

Otherwise `run_reco.sh` will not know that the run belongs to that batch.

---

# 4. Automatic reconstruction of a complete batch

If all the required calibration and data files are available, the easiest way to reconstruct a complete batch is:

```bash
./run_reco.sh BATCH
```

For example:

```bash
./run_reco.sh 3
```

The script reads `run_list.dat`, finds all runs belonging to batch 3 and runs the reconstruction automatically.

Before running it, make sure that:

1. `run_list.dat` contains the correct runs.
2. The corresponding raw data files are present in `data/`.
3. A calibration file exists in `calib/` for every run.
4. The `hadd` executable used inside `run_reco.sh` points to your ROOT installation.

## Important: `hadd`

On some systems I had to explicitly specify the full path to `hadd`.

Find yours with:

```bash
which hadd
```

For example, on macOS with Homebrew it may look like:

```text
/opt/homebrew/Cellar/root/.../bin/hadd
```

Edit the corresponding line in `run_reco.sh` if necessary.

Do not blindly copy my ROOT path. It is almost certainly different on another computer.

---

# 5. Manual reconstruction

For a new run, it is strongly recommended to perform the reconstruction manually the first time.

This makes it easier to check that the calibration, preprocessing and event building are working correctly.

The reconstruction consists of four steps.

---

## Step 1: p.e. calibration

The calibration is performed using:

```text
calib_phe.C
```

on a corresponding:

```text
phe_RUN.root
```

file.

For example, for run `260824`, use the corresponding p.e. calibration data.

The script produces:

```text
calib/phe_260824.dat
```

containing the calibration constants for each channel.

It also produces a PDF in:

```text
pdf/
```

containing the calibration fits.

### IMPORTANT: CHECK THE FITS

Always inspect the calibration PDF.

Do not assume that the calibration worked just because the script finished without errors.

Check that:

- every active channel has a sensible peak;
- the fit follows the single-p.e. peak correctly;
- there are no obviously failed fits;
- calibration values are reasonable compared with neighbouring channels.

A bad p.e. fit produces a wrong charge calibration and therefore wrong results in all subsequent steps.

---

## Step 2: Data preprocessing

Run:

```text
preprocess.C
```

passing:

1. the run number;
2. the batch number.

For example:

```bash
root -l -q 'preprocess.C(260824,3)'
```

for run `260824`, belonging to batch `3`.

The script produces:

```text
out/260824_preprocessed.root
```

and monitoring plots in:

```text
pdf/
```

### Calibration requirement

The preprocessing code looks for a calibration file with the same run number:

```text
calib/phe_RUN.dat
```

For example:

```text
calib/phe_260824.dat
```

must exist before preprocessing run `260824`.

### If no dedicated calibration was taken

A p.e. calibration was not necessarily taken before every physics run.

In this case, use a calibration from the same detector configuration / crystal batch.

Copy it and rename it using the run number you want to analyse.

For example:

```bash
cp calib/phe_260823.dat calib/phe_260824.dat
```

Only do this if the calibration is appropriate for that run.

---

## Step 3: Event building

Run:

```text
event_building.C
```

passing:

1. the run number;
2. the batch number.

For example:

```bash
root -l -q -b 'event_building.C(260824,3)'
```

### Use batch mode

Always use:

```text
-b
```

when running this script.

The event-building code creates hundreds of canvases. Running it interactively is therefore unnecessary and annoying.

The script produces:

```text
out/260824_evt.root
```

containing the reconstructed event information.

It also produces a PDF in:

```text
pdf/
```

with the monitoring plots.

### Check the output

Before moving to the next step, verify that:

```bash
ls -lh out/260824_evt.root
```

actually returns a valid file.

Do not assume that the file exists merely because ROOT finished running.

---

## Step 4: Merge runs belonging to the same batch

If several runs were taken for the same crystal batch, merge the event files before fitting.

For example:

```bash
hadd -f out/batch3.root \
    out/260824_evt.root \
    out/260825_evt.root \
    out/260826_evt.root
```

Use only runs belonging to the same crystal batch.

Do not merge runs from different batches.

The automatic `run_reco.sh` script performs this step automatically.

---

# 6. Light-yield and attenuation-length fits

The final analysis is performed with:

```text
ana_fit.C
```

The script takes the batch number as input.

For example:

```bash
root -l -q -b 'ana_fit.C(3)'
```

Again, batch mode is recommended because the script produces a large number of plots.

The output PDF contains:

- charge fits for the individual crystal positions;
- light-yield estimates;
- light-yield trends along each crystal;
- summary plots;
- attenuation-related information.

### IMPORTANT: CHECK THE FITS

The analysis is automated, but automatic fits can fail.

Always inspect the output PDF before using the fitted values.

In particular, check for:

- fits that clearly do not follow the charge peak;
- unrealistically large parameter uncertainties;
- parameters stuck at their allowed limits;
- anomalous points in the LY-vs-position plots;
- channels with very low statistics.

A ROOT fit returning successfully does not automatically mean that the result is physically meaningful.

---

# 7. Recommended workflow for a new batch

When analysing a new batch for the first time:

```text
1. Copy the raw data
        ↓
2. Rename the files
        ↓
3. Update run_list.dat
        ↓
4. Run / check the p.e. calibration
        ↓
5. Run preprocess.C manually
        ↓
6. Check monitoring plots
        ↓
7. Run event_building.C manually
        ↓
8. Check the output
        ↓
9. Once everything looks correct, use run_reco.sh
        ↓
10. Run ana_fit.C
        ↓
11. Check all final fits
```

Do not immediately run the entire automatic chain on a completely new detector configuration without checking at least one run manually first.

---

# 8. Quick example

For batch 3 and run 260824:

### Copy the data

```bash
scp streamdaq@193.206.147.141:/home/streamdaq/Documents/BGO_test/DAQ/260824/RAW/SDataR_260824.root ./data/
```

### Rename

```bash
mv data/SDataR_260824.root data/260824.root
```

### Check calibration

Make sure that:

```text
calib/phe_260824.dat
```

exists.

### Preprocess

```bash
root -l -q 'preprocess.C(260824,3)'
```

### Event building

```bash
root -l -q -b 'event_building.C(260824,3)'
```

### Check output

```bash
ls -lh out/260824_evt.root
```

### Analyse the full batch

Once all runs belonging to batch 3 are reconstructed:

```bash
root -l -q -b 'ana_fit.C(3)'
```

---

# 9. Common problems

## `run_reco.sh` does not find a run

Check:

```text
run_list.dat
```

and make sure the run is listed under the correct batch.

---

## `preprocess.C` cannot find the calibration

Check that:

```text
calib/phe_RUN.dat
```

exists.

For example:

```bash
ls calib/phe_260824.dat
```

---

## `hadd` is not found

Run:

```bash
which hadd
```

and update the `hadd` path in `run_reco.sh`.

---

## ROOT says a file does not exist

Check the filename and folder manually:

```bash
ls data/
ls calib/
ls out/
```

In particular, check that CoMPASS files have been renamed from:

```text
SDataR_RUN.root
```

to:

```text
RUN.root
```

---

## The code runs but the result looks wrong

Do not trust the output automatically.

Check:

1. the p.e. calibration fits;
2. the preprocessing monitoring plots;
3. the event-building plots;
4. the final LY fits;
5. that the correct batch number was used;
6. that the correct calibration file was used.

If one of the intermediate steps is wrong, rerunning the final fit will not magically fix it.
