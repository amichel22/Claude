# Script: Lin, Lu & Fan 2025 journal club

About 16 minutes when read at a normal pace. [click] means click once for the next item on the same slide. Say "BOLDsωimsuite" as "bold swim suite".

## Slide 1 · Title

Hi. Today I'm presenting the technical review by Chieh-Te Lin, Hanzhang Lu and Audrey Fan, published in Magnetic Resonance in Medicine. It covers how MR fingerprinting is being used to measure brain physiology.

I picked it because it lines up with my rotation project. Since you know this paper better than anyone, please stop me if I get something wrong.

## Slide 2 · Standard methods each cover one part of the picture

The paper starts with a practical problem. Each standard MRI method for brain physiology covers only one part of the picture.

DSC gives blood volume, blood flow and transit time, but it needs a gadolinium bolus, which limits it in patients with kidney problems. ASL measures blood flow without contrast, but the signal is small, so SNR and resolution are low. DCE measures leakage across the blood-brain barrier, again with gadolinium.

So to get a full picture, you stack several scans. That means a longer exam, possible misregistration between scans, and data that are hard to combine.

[click]

Physiological MRF tries to get several of these parameters from one acquisition, matched against one model.

## Slide 3 · MRF lets the signal vary on purpose

Briefly, since this part is familiar: fingerprinting comes from Ma and colleagues in 2013. Instead of fixed settings, MRF changes the flip angle and the repetition time at every time point.

In the animation, the top two strips are the flip angle and the repetition time. Below them are simulated signals for white matter, gray matter and CSF. Each tissue has its own T1 and T2, so each one traces a different curve. That curve is the fingerprint.

To read it, you simulate the same train for many T1 and T2 values and look for the closest curve.

## Slide 4 · Reading a fingerprint: dictionary matching

Here's how that search works. The animation uses a small vascular dictionary I simulated, the same kind of signal we'll see later for MRvF.

The dark dots are a simulated voxel signal with a little noise added, standing in for a measurement. The gray curve is the entry being checked, and the red curve is the best match so far. For each entry, you normalize both curves and take the inner product. That gives a similarity score, and the highest score wins.

In the grid, each square is one combination of oxygen saturation and vessel radius, and darker means more similar. When the search finishes, the best entry is circled, and its parameters become that voxel's values. Here it lands on the true values. So one search gives you all the parameters at once.

Notice that the scores are all above 0.99. Neighboring entries look very alike, which is why noise matters so much.

## Slide 5 · Physiological MRF puts physiology into the dictionary

This is Figure 1 of the paper, and it frames the whole review. The matching idea is the same as standard MRF, but two things change.

First, the acquisition is built to be sensitive to physiology, like blood flow, oxygenation or water exchange, not just T1 and T2.

Second, the dictionary comes from a biophysical model. Here the model is a virtual voxel of vessels, over the ranges of oxygen saturation, blood volume and radius shown at the bottom left.

Each voxel's signal is matched to that dictionary, and you get maps of SO2, CBV and vessel radius.

## Slide 6 · The review sorts the field into three families

The paper sorts the field into three families, and the rest of it follows this order.

MRF-ASL randomizes the labeling and fits a kinetic model to get CBF and bolus arrival time.

MR vascular fingerprinting, or MRvF, gets blood volume, mean vessel radius and SO2 from gradient and spin echo signals, using a Monte Carlo model of vessels. Repeated during a vascular challenge, it also gives cerebrovascular reactivity.

The third, MRF-VFD, for vascular fluid dynamics, targets vessel function and permeability, like Ktrans and water exchange across the barrier.

Many of these methods also give supporting maps like T1, T2 and B1.

## Slide 7 · MRF-ASL: randomize the labeling

Starting with MRF-ASL. This figure is adapted from Su and colleagues.

In panel A, label and control blocks come in a pseudo-random order. In panel B, the labeling duration gets shorter across the scan. Together, one scan samples hundreds of effective delays, instead of the handful you get from multi-delay ASL.

Panel C shows the models. The one-compartment model fits B1, CBF, tissue T1 and tissue arrival time. The two-compartment model adds arterial arrival time, blood volume and travel time.

Wright and colleagues mapped CBF and arrival times in 2.5 minutes, and their CBF agreed with pseudo-continuous ASL, or pCASL.

The bar chart is test-retest data from a later Su paper, about 20 days apart on average. B1 and T1 are very repeatable, above 0.96. Arrival time and CBF are lower, at 0.86 and 0.75. As the paper notes, the parameters that rely on the kinetic model and the small ASL signal are the least repeatable.

## Slide 8 · MRvF: the GESFIDE echo train

Now MR vascular fingerprinting, the part closest to my project. This is the GESFIDE sequence the lab uses, in a figure adapted from Ni and colleagues.

In part A, after the 90-degree pulse, gradient echoes sample the free induction decay. The signal falls with R2 star, which includes dephasing from static fields.

In part B, after the 180-degree pulse at 50 milliseconds, the phase from static fields starts to unwind, and it's refocused at the spin echo at 100 milliseconds. T2 decay keeps going, so the signal still falls, just more slowly. Comparing the slopes in parts A and B is how GESFIDE separates R2 from R2 prime.

After the spin echo, in part C, dephasing builds again.

White matter, gray matter and CSF already follow different curves. The question is how this curve carries vessel information.

## Slide 9 · Why one curve can report SO2, CBV and R

This animation is my attempt to show why. I simulated it for this talk, using the field around a single vessel, with T2 decay turned off.

Deoxyhemoglobin is paramagnetic, so a vessel carrying deoxygenated blood makes a small field around itself. The colors in the background show that field. Each dot is a water spin, colored by its phase.

On the left is a 3-micron vessel, and on the right a 30-micron one, drawn at a tenth the scale. Watch the trace under each one. Both signals drop at first as the spins dephase.

At 50 milliseconds, the 180-degree pulse flips the phases. Around the large vessel, the field barely changes over the distance a spin diffuses, so each spin sees nearly the same field before and after the 180, and the phases line up again. By the spin echo, the signal is back to about 0.92. Around the small vessel, water diffuses through the field, so the 180 only partly undoes the dephasing. The signal ends near 0.75.

So SO2 sets the field strength, CBV sets how much of the voxel feels it, and vessel size sets how much the spin echo recovers.

## Slide 10 · MRvF dictionaries come from virtual voxels

To turn that into numbers, MRvF builds the dictionary from virtual voxels, shown in Figure 3B. You fill a voxel with randomly oriented cylinders and let water diffuse through it. Then you compute the signal for each SO2, CBV and radius.

On the right is the GESFIDE signal from my small dictionary as I change one parameter at a time, again without T2 decay. The gray curve is where it starts. Lower SO2 makes the signal decay faster. Changing the radius shifts the balance between the gradient echo and spin echo parts. More blood volume also speeds up the decay.

Christen's 2014 paper built this kind of dictionary, with about 53,000 entries, and used it with scans before and after ferumoxytol. Newer work needs no contrast: SAGE is fast enough to track reactivity, and Coudert used bSSFP.

In rodent tumors, dictionaries built from realistic 3D vessel networks also matched histology better than cylinders.

## Slide 11 · MRF-VFD: exchange across the blood-brain barrier

The third family, MRF-VFD, looks at exchange across the blood-brain barrier. Panel A, from Venugopal and colleagues, is a hybrid sequence: a gradient-echo readout, a refocusing pulse, then a spin-echo readout. Panel B is a model from Thomson and colleagues: a blood compartment and a tissue compartment that exchange water. The key parameter is tau b, how long water stays in a vessel before it crosses the barrier.

On the right are Thomson's maps, without contrast, in three healthy volunteers. Gray matter has higher CBV than white matter, and tau b and blood T1 look similar across regions.

With contrast, MacAskill's DCE-MRF cut the variability of T1 over time two- to three-fold compared with standard DCE. That was in mice with breast tumors, so it's a tumor perfusion result rather than a brain one.

Acquisition choices also involve trade-offs. Van Dorth's simulations suggested that skipping the preload dose helps permeability estimates but can make blood volume less accurate.

## Slide 12 · Clinical use: stroke and Moyamoya disease

Moving to clinical results. This is MRF-ASL in 34 ischemic stroke patients, from Hongli Fan and colleagues in Hanzhang Lu's group. Each row is one patient, with stroke severity increasing going down. Lesions are outlined in green and normal tissue in blue.

In the lesions, one-compartment CBF was about 15 versus about 29 in normal tissue, and arrival time was about 1.5 seconds versus 0.9. These values tracked stroke severity, and MRF-ASL classified lesions better than pCASL.

In Moyamoya disease, Su and colleagues found longer arrival times behind the narrowed arteries. CBF stayed normal there and went up in the posterior circulation, which the paper reads as collateral flow.

## Slide 13 · Clinical use: tumors and vessel remodeling

For tumors, this figure compares MRvF with vessel size imaging in glioma patients, from Venugopal and colleagues. The radius and blood volume maps look similar, and Venugopal reported MRvF as more robust to noise, with higher resolution.

In rat glioma models, Lemasson and colleagues found that 9L tumors had SO2 around 30 percent, compared with 45 to 50 percent in the C6 and F98 models. Conventional quantitative BOLD, or qBOLD, overestimated it at 40 to 50 percent.

And in atherosclerotic mice, Pouliot found larger vessels, more blood volume and lower oxygen saturation than in wild-type mice. So MRvF picks up vessel remodeling, not only tumors.

## Slide 14 · Machine learning replaces the dictionary search

Next is machine learning, the part I'm most interested in for the rotation.

Dictionary matching compares every voxel with every entry, so the cost grows with the dictionary. A neural network learns the mapping from signal to parameters once, and after that it's a single pass per voxel.

This figure is also from Hongli Fan. The network rows, labeled ANN, are cleaner than dictionary matching, DM, and the green arrows mark artifacts the network removed. Total processing went from 3 hours 12 minutes to 3.6 seconds. DeepMARS reported under half a millisecond per voxel, versus more than 4 seconds.

For MRvF, Lin and colleagues trained a network on simulated GESFIDE dictionaries and got better accuracy than matching for CBV, radius and SO2.

The caveat is that most of these networks learn from simulations. This one is an exception, trained on experimental data. If a voxel looks unlike anything in training, a network can still return a plausible number, with no warning.

## Slide 15 · What still limits physiological MRF

There are still limits, and the paper is clear about them. I'll go through four.

One, field inhomogeneity. Errors in B0 and B1 bias the fit. You can map them and put them in the dictionary, use sequences that are less sensitive to B0, or correct with deep learning.

Two, SNR. Wheeler's simulations suggest MRvF needs an SNR of about 20 or more to match well. In ASL, the label is only about one percent of the signal.

Three, model assumptions. Dictionaries are usually built for healthy tissue and simplified models. In disease, the real signal may not look like any entry, and the match gets biased.

Four, validation. The paper suggests static phantoms, then flow phantoms, then comparison with PET in vivo, plus shared protocols across sites.

Looking ahead, it points to brain water exchange and glymphatic flow, and to simultaneous PET/MR.

## Slide 16 · My take, and questions for discussion

Here's my take. What I found most useful is that the paper puts ASL, vascular and fluid-dynamics MRF into one framework: acquisition, dictionary, matching. Table 2 works as a lookup when you want a study's dictionary ranges or matching method. And it's direct about the limits.

[click]

I also have three questions I'd like to discuss.

First, besides PET, what does the lab use as a reference for vessel radius and SO2 in humans?

Second, the dictionary depends on how it was simulated. A simulation toolbox paper, BOLD swim suite, found that methods differ in radius error. Should dictionaries be benchmarked and versioned the way sequences are?

Third, a network inherits its dictionary's assumptions. Besides physics-informed training, how would we flag voxels outside its range, like a tumor or an infarct?

## Slide 17 · How this connects to my rotation

Last, how this connects to my rotation. The plan follows the MRvF and machine learning sections of the paper.

BOLD swim suite simulates the vessel geometry, the field, diffusion and the signal. I'll use it to build a GESFIDE dictionary with the lab's 38-echo timing. Then I'll compare dictionary matching with a neural network, and apply both to the lab's data.

The questions are simple. Is the network as accurate? How much faster is it? And how does it behave at low SNR, or outside its training range?

My first steps are to run Jack's notebook, reproduce the reference simulation, and rebuild my test dictionary in BOLD swim suite.

## Slide 18 · Thank you

That's everything. Thank you, and I'm happy to take questions.

## If asked

**Why does the signal still fall after the 180 in GESFIDE?** Between the 180 and the spin echo, static dephasing is being reversed, but T2 decay continues, so the signal falls at R2 minus R2 prime instead of R2 plus R2 prime. The slope difference between parts A and B gives R2 prime. If the drop across the 180 itself looks large, that is likely imperfect refocusing, for example a slice-profile mismatch between the 90 and the 180.

**Who introduced GESFIDE?** Ma and Wehrli, in 1996. The figure in the review is adapted from Ni and colleagues' 2015 comparison of R2 prime methods at 3T.

**Why GESFIDE and not SAGE?** GESFIDE samples many echoes before and after one 180, so you see the free induction decay, the refocusing and the spin echo in detail. The paper notes it separates T2 and T2 star more directly and stably than bSSFP, without banding, at the cost of lower SNR and less susceptibility contrast. SAGE uses fewer echoes but is fast enough for dynamic maps, like Wheeler's reactivity work.

**Why did Christen use the ratio of pre- and post-contrast signals?** The ratio cancels anything that scales both scans the same way: proton density, coil sensitivity, T2 decay, and macroscopic B0 and B1 effects. What's left is mainly the extra dephasing the contrast agent adds, which depends on CBV and radius. SO2 doesn't cancel exactly, because the dephasing isn't linear in susceptibility and the dictionary simulates both states with deoxyhemoglobin included, so the same match also gives SO2.

**What is tau b?** The average time a water molecule stays inside a vessel before it crosses the blood-brain barrier. Higher water permeability gives a shorter tau b. Its inverse is the water exchange rate, roughly the permeability-surface area product for water divided by capillary blood volume. It's about water, so it can change even when the barrier doesn't leak contrast agent.

**Why the inner product and not least squares?** Normalizing both curves removes the unknown overall scale, like proton density and coil gain, so you only compare shape. For unit-length vectors, the squared difference is 2 minus 2 times the inner product, so the highest inner product is the smallest error. It's also one fast matrix multiply.

**Why is a network faster, and what can it get wrong?** Matching compares each voxel with every entry, so the time grows with the dictionary, and the answer can only land on a grid point. A trained network does a fixed amount of work per voxel and gives continuous values. The risk is input unlike its training data: it can return a plausible but wrong value with no flag, so it has to be tested on other scanners, patient groups, known artifacts and extreme values.

**What does the two-compartment ASL model add?** It separates labeled blood still passing through arteries from blood that has reached tissue. That adds arterial arrival time, arterial blood volume and travel time, which matters when arrival is delayed, as in stroke or Moyamoya.

**Why does MRvF need an SNR of about 20?** Neighboring dictionary entries differ by small amounts. In Wheeler's simulations, matching errors for SO2, CBV and radius dropped clearly once SNR reached about 20.

**Why is CBF less repeatable than T1 in MRF-ASL?** CBF and arrival time come from fitting a kinetic model to a label signal that is about one percent of the total. T1 and B1 come from the much larger tissue signal.

**The Moyamoya numbers.** Table 3 lists arrival time as 1900 versus 1100 milliseconds and posterior-circulation CBF as 50 versus 35. The text gives about 1470 versus 964 for arrival time and 40 versus 31 for posterior CBF. I wasn't sure which analysis each set came from, so I kept the qualitative result on the slide.

**How did you make the animations?** All four are my own simulations in Python, for illustration only. The fingerprint curves come from an extended phase graph simulation of an inversion-recovery FISP train with 500 time points. The spin animation is a single vessel perpendicular to B0 at 3 tesla, SO2 60 percent, hematocrit 0.4, diffusion 1 square micron per millisecond, 40,000 spins, T2 off. The dictionary is a 2D Monte Carlo model with randomly placed cylinders, each with its own random B0 direction, using the lab's 38-echo timing, over SO2 30 to 95 percent, radius 2 to 30 microns and CBV 1 to 6 percent. The code is in my repo.

**How does BOLD swim suite fit in?** It is an open Python toolbox for forward modeling the BOLD signal: vessel geometry, field offsets, diffusion and the pulse sequence. For MRvF it would generate the dictionary. Its paper shows that simulation choices change radius estimates, which is why I asked about benchmarking dictionaries.
