# Toolbox talk: read-aloud script

One section per slide, in deck order. The same text is in each slide's speaker notes.

## Slide 1 · BOLDsωimsuite

Hi everyone. Today I want to walk through the toolbox Jack keeps mentioning, BOLD swim suite.

I'll cover four things: what it does and why the lab needs it, the two papers it comes from, what I ran and played with, and how I'd turn it into my rotation project.

I'm still new to MRI physics, so please stop me if I get something wrong.

## Slide 2 · Why the lab needs a simulator

First, why we need a simulator at all.

With vascular fingerprinting, the scanner gives us a decay curve for each voxel, 38 echoes from GESFIDE. On its own, that curve doesn't tell us blood volume or vessel size.

So we build a dictionary: simulated curves for every combination of oxygenation, blood volume and vessel radius. Then we match each measured curve to the closest simulated one and read off the parameters. That gives us the maps.

The key point is on the bottom. The maps can only be as right as the simulated curves. If the simulation has a bias, the maps inherit it in a systematic way, and the matching step has no way of flagging it. The curve still matches well; it's just the wrong curve.

BOLD swim suite is the tool that makes those simulated curves.

## Slide 3 · What one simulation does

Here's what a single simulation looks like, from the toolbox paper.

On the left is a virtual voxel: a 2D slice with vessels crossing it, like a capillary bed. Each vessel gets its own random field direction, which stands in for vessels pointing every which way.

In the middle is the magnetic field around them. Deoxygenated blood is slightly magnetic, so each vessel bends the field nearby. That's the red and blue pattern.

Then the simulation lets water molecules wander through that field. Each one picks up a slightly different phase, they fall out of sync, and the signal drops. On the right you can see the spin-echo signal: it decays, the refocusing pulse brings part of it back, and then it decays again. They also split it into the tissue part and the blood part.

For a dictionary, you just repeat this for every combination of oxygenation, blood volume and radius.

## Slide 4 · Every simulation is four choices

This is the figure that made the toolbox click for me. Every simulation is four choices, left to right.

Geometry: is the voxel 2D or 3D, and are the vessels straight cylinders, spheres, or a realistic vessel network.

The field: do you use the exact formula for a cylinder, or a Fourier transform on a grid, which works for any shape but is approximate.

Diffusion: do you track thousands of random water molecules, which is Monte Carlo, or blur the magnetization on a grid, which is deterministic diffusion.

And the signal: which pulse sequence you play, for us the 38-echo GESFIDE.

The method name just spells out the choices. 3D-ANA-MC means 3D, analytical field, Monte Carlo. MrVox, the older Matlab tool, is one fixed combination: 2D, FFT, grid diffusion, plus a shortcut called 3B0 that I'll come back to, because it turns out to matter.

## Slide 5 · Paper 1: the toolbox, and what each choice costs

Now the two papers Jack keeps pointing me to. The first one is the toolbox paper, from Jacob Chaussé, Avery Berman and Jean Chen.

Their point is that BOLD simulations used to live in separate codes that mostly weren't shared. They put the common methods into one Python package, free on GitHub, with tutorial notebooks.

They also measured what each choice costs. 2D is faster, 3D is more realistic. A grid voxel runs faster per step, but in 3D it takes about a hundred times more memory. The FFT field works for any shape but isn't exact. And grid diffusion gives a clean, noise-free answer, but gets very memory hungry for big vessels. The figure at the bottom shows how big these costs get: a 3D grid voxel took about an hour to build and fourteen and a half gigabytes of memory, versus a tenth of a gigabyte for the continuous one.

Their general advice: use the exact field when you can, and keep enough grid points across each vessel.

## Slide 6 · Paper 2: does the recipe change the answer?

The second paper uses the toolbox to ask the question that matters for us: does the recipe change the answer?

They took eight approaches people have published and ran all of them with the same settings: 3 tesla, 2 percent blood volume, fifteen vessel sizes from 1 to 60 microns, and ten random voxels.

The figure shows the geometries. From left: a single vessel in 2D, many vessels in 2D, either each with its own field direction or with a three-direction average, random cylinders in 3D, and a simplified vessel network with branching capillaries.

The reference is 3D cylinders with the exact field and Monte Carlo diffusion.

They compared everything in three ways, and the one closest to what we do is a mini fingerprinting test: use each approach as the dictionary and see if you get the right vessel radius back.

## Slide 7 · Two shortcuts got the vessel size wrong

Here's the fingerprinting result. Each small panel uses one approach as the dictionary. The horizontal axis is the true radius, the vertical axis is the radius the matching picked, and the red squares are the estimates. On the black line means it got it right.

Most approaches are on the line. The 2D version with a random field direction per vessel got all fifteen right, and so did grid diffusion. The FFT approaches missed two each.

Two shortcuts did badly. The single-vessel approach missed seven, with a mean error of 13 percent. And 3B0, top right, missed twelve of fifteen, with a mean error of 54 percent and a worst case of 233 percent.

So the recipe really can change the vessel size you report.

## Slide 8 · Why 3B0 fails: the blood stops dephasing

Why does 3B0 do so badly? I made these two field maps to see it for myself. Same vessels, two ways of computing the field.

On the left, each vessel gets its own field direction, which is what happens in a real brain, since vessels point every which way. The patterns point in different directions, and the blood inside each vessel sits at a different field.

On the right is 3B0. It's a 2D shortcut that averages the field over three directions to mimic vessels pointing every which way. Every pattern now points the same way, and the blood in every vessel sits at the same offset, close to zero.

If all the blood sits at the same offset, there's no spread between vessels to dephase, so the blood signal hardly decays. In my runs it decays at 2 to 3 per second, versus about 42 in the 3D reference. Berman shows the same thing.

That sounds small, because blood is only a few percent of the voxel. But blood normally decays so fast that it carries at least a quarter of the gradient-echo decay. Berman found about 25 percent at 2 percent blood volume, and in my runs it's over half for capillary-sized vessels. So the whole signal is off, and so is anything you fit from it.

## Slide 9 · Why this matters for our dictionaries

So why does this matter for us?

MrVox, the Matlab toolbox, uses exactly that 3B0 shortcut, together with an FFT field and grid diffusion. It made the dictionary for the first vascular fingerprinting paper, Christen 2014, and several after it. As far as I understand, our current dictionaries come from MrVox too. Jack, correct me if that's wrong.

And it's not only Berman saying this. People who work on MrVox compared 3B0 with a dictionary from realistic vessel networks in rats, and 3B0 overestimated vessel radius by 14 to 19 percent.

Berman's conclusion is that the 3D reference is a good default, and the 2D version where each vessel gets its own field direction reproduced it almost exactly. They say 3B0 should be used with caution.

So what's missing? Berman only tested radius, but with a known answer. Delphin's group did use GESFIDE and estimated all three parameters, but in rats, where you don't know the true values. A known-answer test of 3B0 on our 38-echo GESFIDE, for blood volume, radius and oxygenation together, hasn't been done. That's the gap I think I can work on.

## Slide 10 · What I ran: the install, the examples and Jack's notebook

Now what I actually ran.

The toolbox only works with Python 3.10, so I made a separate environment for it. Then I ran Jack's install check, the example script that comes with the toolbox, and Jack's GESFIDE notebook from start to end.

The runs are seeded, so the numbers should repeat exactly, and mine matched Jack's. A few examples on the right: the field next to a vessel matches the textbook formula. Water spreads 39.9 square microns in 20 milliseconds, against 40 from theory. The GESFIDE fit on the reference voxel gives R2 of 14.8 and R2 prime of 1.81 per second.

And the last row is a preview of the problem. In Jack's head-to-head on our sequence, one voxel each, the blood part of the signal from the MrVox-like 3B0 voxel has an NRMSD of 94 percent against the 3D reference. That's Berman's error measure, the RMS difference relative to the reference's spread. The 2D version with its own field direction per vessel is at about 11 percent.

## Slide 11 · Playing with it: what each knob does to our signal

Then I started playing with it on our own sequence: the 38 GESFIDE echoes, with the 180 pulse at 50 milliseconds and the spin echo at 100.

On the left I change oxygenation. Less oxygen means more magnetic blood, so the signal decays faster. That's the BOLD effect.

On the right I change vessel size, and the shape changes. Big vessels behave like a static field, so the 180 pulse brings a lot of signal back at the spin echo. Small vessels don't refocus, because water diffuses past them too fast. That difference in shape is what lets fingerprinting tell vessel sizes apart.

A nice detail: I only ran one simulation per vessel size and voxel. Oxygenation just scales the field, so one run gives every SO2 for free. Jack's notebook shows that trick, and it makes the dictionary much cheaper.

If someone asks why most curves peak a bit before 100 milliseconds: diffusion causes some loss that the 180 can't undo, and that loss keeps growing, so the maximum shifts earlier. The 30 micron curve, where diffusion matters least, peaks closest to 100.

## Slide 12 · Berman's result holds on our sequence

Next I wanted to see whether Berman's result holds on our sequence. These are Boxerman plots, the classic check: how fast the signal decays as a function of vessel size.

On the left is the gradient-echo part, at about 30 milliseconds. The 3D reference, in black, and the 2D version with its own direction per vessel, in blue, sit on top of each other. 3B0, in red, sits well below: about three quarters of the reference for large vessels, and much less for small ones.

On the right is the spin echo, at 100 milliseconds. It peaks around capillary size, which is the classic result, and here the three are much closer, with 3B0 a bit low for small vessels.

So the problem is in the reversible part of the decay, R2 prime, which weights every GESFIDE echo except the ones near the 100 millisecond spin echo. Same pattern as Berman's Figure 8C, now with our timing.

## Slide 13 · First try: which dictionary gets the radius right?

This is the part I'm most excited about: a first, small version of the project.

I treated 3D reference simulations as if they were measured data, and matched them against three dictionaries, built three ways. Each panel is one dictionary. Horizontal is the true radius, vertical is what the matching picked. On the gray line means right.

On the left, a 3D dictionary built from new random voxels: 10 of 12 right, and the misses are just the next grid value. In the middle, the 2D recipe with its own direction per vessel: 8 of 12, again only one step off. On the right, the MrVox-style 3B0 dictionary: 3 of 12. It calls a 1 micron vessel 15 microns, and everything small comes out around 7.

That's Berman's finding again, but on our sequence, with oxygenation estimated at the same time. Adding noise at an SNR of 100 didn't change the picture.

To be clear, this is a toy: fixed blood volume, no T2, and only two voxels per entry. It also uses the total signal, blood plus tissue. Most of the 3B0 error lives in the blood signal; in my runs the tissue-only decay rates were only 3 to 16 percent low. So whether our dictionary includes blood matters a lot, and that's one of my questions at the end.

## Slide 14 · The same matches also shift oxygenation

Because the matching picks radius and oxygenation together, the wrong radius drags oxygenation with it.

The true SO2 here is 60 percent. The two good recipes stay within about 7 points for vessels of 2 microns and up, and mostly within 5.

3B0, in red, reads small vessels as up to 20 points too oxygenated, and large vessels 8 to 12 points too low. That's a big error for a number we'd want to compare between patients.

One honest caveat: even the 3D control comes out a few points low. I think that's noise from using only two random voxels per entry, since Berman used ten. That's one of the first things I'd check once I can run bigger jobs on the cluster.

## Slide 15 · Scaling up: why the cluster

Everything so far ran on a laptop. A real dictionary is bigger, and that's where the cluster comes in.

One 3D simulation takes about 35 seconds on one core at 3 percent blood volume; more vessels take longer. A first real dictionary, say 25 radii, 15 blood volumes and 3 random voxels each, is about 1,100 runs. That's 10 to 12 hours on one core if every entry cost the same, and more once the high blood volumes are in. But each run is independent, so as 30 jobs side by side on the cluster it takes about an hour.

Two things make this manageable. Oxygenation is free, because one run can be rescaled to any SO2. And each run uses one core, so it maps onto array jobs the same way the MrVox jobs are already split.

My next step there is to get my account, test one job interactively, and then run a small array with Jack's template.

## Slide 16 · How I would turn this into my rotation

So here's how I'd turn this into my rotation.

The question is in two halves. First, does the way we simulate the dictionary change the blood volume, radius and oxygenation we report from GESFIDE? Second, once the dictionary is right, can a neural network replace the matching step? The second half is what's in my rotation agreement. The first half is what I'd add, because the network can only be as good as the dictionary it learns from.

Part one, over the next two weeks: build a GESFIDE dictionary the recommended way. Small on my laptop first, then full size on the cluster. It's done when it reproduces Jack's reference.

Part two, about three weeks: build a 3B0 dictionary with the same settings and measure the difference in all three parameters, with noise, and then on one lab dataset.

Part three, if there's time: a small network against matching.

The last week is the summary, the lab-meeting talk, and a clean repo. Even if part three doesn't happen, parts one and two are a complete result.

## Slide 17 · Where I would like your input

I'd really like your input on a few things.

First, is "does the simulation recipe change our estimates" a good first question, before the machine-learning part? Or would you rather I go straight to the network?

Second, what would help the lab most: a new dictionary you'd actually adopt, or a measured bias for the current one?

Third, which dataset should I use, and was it acquired with or without contrast agent? That changes what has to go into the dictionary.

And fourth, should the new dictionary include the blood signal, or tissue only? Given what Berman found about the blood signal, that choice matters.

## Slide 18 · Thank you

That's it. Thank you, and I'm happy to take questions.
