# Toolbox talk: read-aloud script

One section per slide, in deck order. The same text is in each slide's speaker notes.

## Slide 1 · BOLDsωimsuite

Hi everyone. Today I want to walk through the toolbox Jack keeps mentioning, BOLD swim suite.

I'll cover four things: what it does and why the lab needs it, the two papers it comes from, what I ran and played with, and how I'd turn it into my rotation project.

I'm still new to MRI physics, so please stop me if I get something wrong.

## Slide 2 · Follow-up from last week

Before I go on, I want to answer the three questions from last week. I spent this week making sure I understand them.

First, you asked whether MrVox computes the field analytically or with an FFT. It's the FFT. MrVox works on a 2D grid, gets the field with a Fourier method, moves the water with grid diffusion, and averages the field over three B0 directions, the 3B0 shortcut. Christen's 2014 paper describes the Fourier field and the grid diffusion, and Chaussé and Berman both say MrVox also uses the 3B0 average. Analytical or FFT only describes how the field is computed; MrVox is a numerical simulation either way. One subtlety: Berman's 3B0 test used the exact field and Monte Carlo, so the errors on the 3B0 panel come from the averaging, not from the FFT. The FFT recipes only missed two radii each.

Second, the y axis on Berman's Figure 4, and why it makes sense to compare these models at all. I made two new slides for that, 8 and 9. Since we covered slides 3 to 7 last week, I'll go straight there after this.

Third, your comment about the purple, the orange and the red. I wasn't sure which figure you meant, so I went through every figure. On slide 4 the coloured dots only label different vessels; the caption says the colours encode different perturbers. Two figures use those colours to mean something. In Berman's Figure 6, purple is the blood signal, orange the tissue signal and green the total. In Figure 8, red-orange and yellow-orange are the two full 2D recipes and purple is the single-vessel shortcut. I'll show Figure 6 on slide 11, because it's about blood versus tissue, and mention Figure 8 there. I also tested the idea with my own simulations: with the tissue signal only, the 3B0 dictionary gets 10 of 12 radii right instead of 3. Could you tell me which figure you meant, so I know I understood your point?

## Slide 3 · Why the lab needs a simulator

First, why we need a simulator at all.

With vascular fingerprinting, the scanner gives us a decay curve for each voxel, 38 echoes from GESFIDE. On its own, that curve doesn't tell us blood volume or vessel size.

So we build a dictionary: simulated curves for every combination of oxygenation, blood volume and vessel radius. Then we match each measured curve to the closest simulated one and read off the parameters. That gives us the maps.

The key point is on the bottom. The maps can only be as right as the simulated curves. If the simulation has a bias, the maps inherit it in a systematic way, and the matching step has no way of flagging it. The curve still matches well; it's just the wrong curve.

BOLD swim suite is the tool that makes those simulated curves.

## Slide 4 · What one simulation does

Here's what a single simulation looks like, from the toolbox paper.

On the left is a virtual voxel: a 2D slice with vessels crossing it, like a capillary bed. Each vessel gets its own random field direction, which stands in for vessels pointing every which way.

In the middle is the magnetic field around them. Deoxygenated blood is slightly magnetic, so each vessel bends the field nearby. That's the red and blue pattern.

Then the simulation lets water molecules wander through that field. Each one picks up a slightly different phase, they fall out of sync, and the signal drops. On the right you can see the spin-echo signal: it decays, the refocusing pulse brings part of it back, and then it decays again. They also split it into the tissue part and the blood part.

For a dictionary, you just repeat this for every combination of oxygenation, blood volume and radius.

## Slide 5 · Every simulation is four choices

This is the figure that made the toolbox click for me. Every simulation is four choices, left to right.

Geometry: is the voxel 2D or 3D, and are the vessels straight cylinders, spheres, or a realistic vessel network.

The field: do you use the exact formula for a cylinder, which is the analytical option, or a Fourier transform on a grid, the FFT option, which works for any shape but is approximate.

Diffusion: do you track thousands of random water molecules, which is Monte Carlo, or blur the magnetization on a grid, which is deterministic diffusion.

And the signal: which pulse sequence you play, for us the 38-echo GESFIDE.

The method name just spells out the choices. 3D-ANA-MC means 3D, analytical field, Monte Carlo. That's Berman's reference.

And to answer your question from last week directly: MrVox is the FFT kind. In this naming it's 2D-GRD-CYL-FFT-DD: a 2D grid, cylinders, the field by FFT, and grid diffusion, plus the 3B0 shortcut, which averages the field over three B0 directions. That last part is the one that turns out to matter.

## Slide 6 · Paper 1: the toolbox, and what each choice costs

Now the two papers Jack keeps pointing me to. The first one is the toolbox paper, from Jacob Chaussé, Avery Berman and Jean Chen.

Their point is that BOLD simulations used to live in separate codes that mostly weren't shared. They put the common methods into one Python package, free on GitHub, with tutorial notebooks.

They also measured what each choice costs. 2D is faster, 3D is more realistic. A grid voxel runs faster per step, but in 3D it takes about a hundred times more memory. The FFT field works for any shape but isn't exact. And grid diffusion gives a clean, noise-free answer, but gets very memory hungry for big vessels. The figure at the bottom shows how big these costs get: a 3D grid voxel took about an hour to build and fourteen and a half gigabytes of memory, versus a tenth of a gigabyte for the continuous one.

Their general advice: use the exact field when you can, and keep enough grid points across each vessel.

## Slide 7 · Paper 2: does the recipe change the answer?

The second paper uses the toolbox to ask the question that matters for us: does the recipe change the answer?

They took eight approaches people have published and ran all of them with the same settings: 3 tesla, 2 percent blood volume, fifteen vessel sizes from 1 to 60 microns, and ten random voxels.

The figure shows the geometries. From left: a single vessel in 2D, many vessels in 2D, either each with its own field direction or with a three-direction average, random cylinders in 3D, and a simplified vessel network with branching capillaries.

The reference is 3D cylinders with the exact field and Monte Carlo diffusion.

They compared everything in three ways, and the one closest to what we do is a mini fingerprinting test: use each approach as the dictionary and see if you get the right vessel radius back.

## Slide 8 · Why compare simulation recipes at all?

You asked why it makes sense to compare these models at all. Here's how I think about it.

In a human scan we can't measure the true vessel radius voxel by voxel. In a simulation we set it, so we know it exactly. So Berman's group made test signals with the 3D reference recipe, where the radius is known, and then asked each of the other recipes to recover it, each using a dictionary built with that recipe. Same physiology and the same sequence; only the recipe changes. So if the reported radius moves, it's the recipe that moved it.

Why should we care? Because a dictionary carries the physics of whatever simulated it, and matching always finds some best entry, even when that physics is off. In Berman's 3B0 panel, several of the wrong picks sit on pale, high-scoring cells, so a good fit didn't mean a right radius. So if two reasonable recipes report different radii for the same signal, the recipe itself is a source of error in every map we make with it.

Why that reference? Random 3D cylinders have an exact formula for the field and no grid artefacts. Berman says Monte Carlo diffusion may be considered the most accurate way to model diffusion. The field has a lot of experience with this model, and many of the later shortcuts were derived from it.

One honest caveat: it's still a model, straight, randomly oriented cylinders. Berman calls it ground truth only for cylinder simulations. So this measures agreement with the best-understood recipe, not the truth in real tissue. But it tells us to check which recipe built our own dictionary.

## Slide 9 · How to read Berman's Figure 4

Here's how to read Berman's Figure 4. I've taken two of its eight panels and numbered the parts.

On the left is a dictionary that works, 2D-ANA-MC.

Number 2, the x axis, is the true radius of the test signal. The test signals were simulated with the 3D reference, so for every column we know the right answer, at least for that model.

Number 1, the y axis, is the radius of each entry in the dictionary. The figure labels it "estimated radius", but the caption says it's the dictionary radius used for comparison. So each column is one test signal compared with every entry in the dictionary.

Number 3, the colour, is how alike the two decay curves are. They used R squared, the same measure Christen 2014 used, after dividing each curve by its average over time, so only the shape counts. The review from your lab, Lin 2025, lists Christen's matching as an inner product; my own matching on slide 17 uses the normalized inner product, which also compares shape only. Near white means almost identical; dark navy means 0.7 or less.

Number 4, the red square, is the best match in each column, the entry with the highest R squared. Its height is the radius the matching reports.

Number 5, the black line, is where the reported radius equals the true radius. Squares on the line are right answers.

Notice that the light band is narrow for the smallest vessels and widens from about 3 to 4 microns up: for bigger vessels, neighbouring sizes give very similar curves. Berman notes that for the FFT recipes at the largest radii, R squared was above 0.99 even between neighbouring sizes, and calls the FFT's small field errors the probable cause of their misses there.

On the right is the 3B0 dictionary. Most squares sit above the line, and every miss is too large: 12 of 15 radii wrong.

One detail: Berman's test signals were a simple spin echo, not our GESFIDE. Running this kind of test on our sequence is what my own simulations add.

## Slide 10 · Two shortcuts got the vessel size wrong

Now the full figure, read the same way. Each panel uses a different recipe to build the dictionary; the test signals are always the same 3D reference signals. The 2D recipes are on top, the 3D ones on the bottom.

The bottom-left panel, 3D-ANA-MC, is the reference matched against its own approach, so a perfect result is expected. It's a check that the matching works, not a finding.

Most approaches are on the line. The 2D version with a random field direction per vessel got all fifteen right, and so did grid diffusion and the gridded 3D version. The two FFT approaches missed two each, mostly at the largest radii, where the curves are hardest to tell apart; Berman calls the FFT's small field errors the probable cause.

Two shortcuts did badly. The single-vessel approach missed seven, with a mean error of 13 percent. And 3B0, top right, missed twelve of fifteen, every one too large, with a mean error of 54 percent and a worst case of 233 percent. With only fifteen sizes in the dictionary, even one step off is a 17 to 50 percent error, which is why these percentages look big.

So the recipe alone can change the vessel size you report. And this is an idealized test: no noise, and only the radius changes. In real data, a bias like this would be harder to notice, not easier.

## Slide 11 · Green, orange, purple: total, tissue and blood

Now the colours. Last week you mentioned the purple. This is Berman's Figure 6, with two of its panels side by side.

The y axis is NRMSD: for each vessel radius, how far a recipe's signal is from the 3D reference over the whole decay, as a percentage of the reference signal's spread over time, its standard deviation. Lower means closer, and it's a log scale.

Each panel shows the same recipe three ways. Green is the total signal. Orange is only the tissue part, the extravascular signal. Purple is only the blood part, the intravascular signal.

On the left, 2D-ANA-MC, green and orange sit together. Purple is higher, around 5 to 10 percent: the blood signal is the hardest part to reproduce. But blood is only 2 percent of the voxel, so that error hardly moves the total. That's why green follows orange.

On the right is 3B0. Purple is about 150 percent off, because the blood in every vessel gets the same field offset, so it barely dephases. Even at 2 percent of the voxel that matters, because in the reference the blood carries about a quarter of the gradient-echo decay. So the green total goes from about 120 percent at the smallest vessels down to about 40, far above the orange tissue curve at 10 to 25 percent. The blood signal is what breaks 3B0.

I think that's what you meant: when green and orange agree, the purple doesn't matter much, and 3B0 is the exception.

There's no red in this figure; the curve that follows orange is green. Berman's Figure 8 is the one with red, orange and purple together. There, the red-orange and yellow-orange markers are the two full 2D recipes, and they sit almost on the black reference line; the purple markers are the single-vessel shortcut, and they don't. Which of the two figures did you mean?

## Slide 12 · Why 3B0 fails: the blood stops dephasing

Why does 3B0 do so badly? I made these two field maps to see it for myself. Same vessels, two ways of computing the field.

On the left, each vessel gets its own field direction, which is what happens in a real brain, since vessels point every which way. The patterns point in different directions, and the blood inside each vessel sits at a different field.

On the right is 3B0. It's a 2D shortcut that averages the field over three directions to mimic vessels pointing every which way. Every pattern now points the same way, and the blood in every vessel sits at the same offset, close to zero.

If all the blood sits at the same offset, there's no spread between vessels to dephase, so the blood signal hardly decays. In my runs it decays at 2 to 3 per second, versus about 42 in the 3D reference. Berman shows the same thing.

That sounds small, because blood is only a few percent of the voxel. But blood normally decays so fast that it carries at least a quarter of the gradient-echo decay. Berman found about 25 percent at 2 percent blood volume, and in my runs it's over half for capillary-sized vessels. So the whole signal is off, and so is anything you fit from it.

## Slide 13 · Why this matters for our dictionaries

So why does this matter for us?

MrVox, the Matlab toolbox, uses exactly that 3B0 shortcut, together with an FFT field and deterministic, grid-based diffusion. It made the dictionary for the first vascular fingerprinting paper, Christen 2014, and several after it, although Christen's paper itself doesn't mention the 3B0 averaging. As far as I understand, our current dictionaries come from MrVox too. Jack, correct me if that's wrong.

To be precise about what Berman tested: their 3B0 version used the exact field and Monte Carlo, so it isolates the averaging step. Separately, the FFT recipes missed only two radii each, and deterministic diffusion alone got all fifteen right. How the three combine inside MrVox itself wasn't tested in these papers. Chaussé mentions an earlier conference abstract from Berman's group that compared MrVox's 3B0 method with other methods and found substantial differences, but I haven't read it.

And it's not only Berman saying this. People who work on MrVox compared 3B0 with a dictionary from realistic vessel networks in rats, and 3B0 overestimated vessel radius by 14 to 19 percent.

Berman's conclusion is that the 3D reference is a good default, and the 2D version where each vessel gets its own field direction reproduced it almost exactly. They say 3B0 should be used with caution.

So what's missing? Berman only tested radius, but with a known answer. Delphin's group did use GESFIDE and estimated all three parameters, but in rats, where you don't know the true values. A known-answer test of 3B0 on our 38-echo GESFIDE, for blood volume, radius and oxygenation together, hasn't been done. That's the gap I think I can work on.

## Slide 14 · What I ran: the install, the examples and Jack's notebook

Now what I actually ran.

The toolbox only works with Python 3.10, so I made a separate environment for it. Then I ran Jack's install check, the example script that comes with the toolbox, and Jack's GESFIDE notebook from start to end.

The runs are seeded, so the numbers should repeat exactly, and mine matched Jack's. A few examples on the right: the field next to a vessel matches the textbook formula. Water spreads 39.9 square microns in 20 milliseconds, against 40 from theory. The GESFIDE fit on the reference voxel gives R2 of 14.8 and R2 prime of 1.81 per second.

And the last row is a preview of the problem. In Jack's head-to-head on our sequence, one voxel each, the blood part of the signal from the MrVox-like 3B0 voxel has an NRMSD of 94 percent against the 3D reference. That's Berman's error measure, the RMS difference relative to the reference's spread. The 2D version with its own field direction per vessel is at about 11 percent.

## Slide 15 · Playing with it: what each knob does to our signal

Then I started playing with it on our own sequence: the 38 GESFIDE echoes, with the 180 pulse at 50 milliseconds and the spin echo at 100.

On the left I change oxygenation. Less oxygen means more magnetic blood, so the signal decays faster. That's the BOLD effect.

On the right I change vessel size, and the shape changes. Big vessels behave like a static field, so the 180 pulse brings a lot of signal back at the spin echo. Small vessels don't refocus, because water diffuses past them too fast. That difference in shape is what lets fingerprinting tell vessel sizes apart.

A nice detail: I only ran one simulation per vessel size and voxel. Oxygenation just scales the field, so one run gives every SO2 for free. Jack's notebook shows that trick, and it makes the dictionary much cheaper.

If someone asks why most curves peak a bit before 100 milliseconds: diffusion causes some loss that the 180 can't undo, and that loss keeps growing, so the maximum shifts earlier. The 30 micron curve, where diffusion matters least, peaks closest to 100.

## Slide 16 · Berman's result holds on our sequence

Next I wanted to see whether Berman's result holds on our sequence. These are Boxerman plots, the classic check: how fast the signal decays as a function of vessel size.

On the left is the gradient-echo part, at about 30 milliseconds. The 3D reference, in black, and the 2D version with its own direction per vessel, in blue, sit on top of each other. 3B0, in red, sits well below: about three quarters of the reference for large vessels, and much less for small ones.

On the right is the spin echo, at 100 milliseconds. It peaks around capillary size, which is the classic result, and here the three are much closer, with 3B0 a bit low for small vessels.

So the problem is in the reversible part of the decay, R2 prime, which weights every GESFIDE echo except the ones near the 100 millisecond spin echo. Same pattern as Berman's Figure 8C, now with our timing.

## Slide 17 · First try: which dictionary gets the radius right?

This is the part I'm most excited about: a first, small version of the project.

I treated 3D reference simulations as if they were measured data, and matched them against three dictionaries, built three ways. Each panel is one dictionary. Horizontal is the true radius, vertical is what the matching picked. On the gray line means right.

On the left, a 3D dictionary built from new random voxels: 10 of 12 right, and the misses are just the next grid value. In the middle, the 2D recipe with its own direction per vessel: 8 of 12, again only one step off. On the right, the MrVox-style 3B0 dictionary: 3 of 12. It calls a 1 micron vessel 15 microns, and everything small comes out around 7.

That's Berman's finding again, but on our sequence, with oxygenation estimated at the same time. Adding noise at an SNR of 100 didn't change the picture.

To be clear, this is a toy: fixed blood volume, no T2, and only two voxels per entry. It also uses the total signal, blood plus tissue. Most of the 3B0 error lives in the blood signal; in my runs the tissue-only decay rates were only 3 to 16 percent low. So whether our dictionary includes blood matters a lot, and that's one of my questions at the end.

## Slide 18 · The same matches also shift oxygenation

Because the matching picks radius and oxygenation together, the wrong radius drags oxygenation with it.

The true SO2 here is 60 percent. The two good recipes stay within about 7 points for vessels of 2 microns and up, and mostly within 5.

3B0, in red, reads small vessels as up to 20 points too oxygenated, and large vessels 8 to 12 points too low. That's a big error for a number we'd want to compare between patients.

One honest caveat: even the 3D control comes out a few points low. I think that's noise from using only two random voxels per entry, since Berman used ten. That's one of the first things I'd check once I can run bigger jobs on the cluster.

## Slide 19 · Without the blood signal, 3B0 mostly recovers

So this week I tested that idea, without running anything new. Each of my 96 simulations had already saved the tissue signal and the blood signal separately, so I redid the matching twice: once on the total signal, as before, and once on the tissue signal only, for the test signals and the dictionaries.

On the left is the total signal, the same result as two slides ago: the 3B0 dictionary gets 3 of 12 radii right.

On the right is tissue only. Now 3B0 gets 10 of 12 right, and the 2D recipe with its own direction per vessel gets 11 of 12. SO2 gets much better too: 3B0 is within 4 points from 3 microns up. The smallest vessels, 1.5 and 2 microns, are still 15 to 20 points off in SO2.

So most of the 3B0 problem lives in the blood signal.

Two honest limits. First, I also took the blood out of the test signals, and we can't do that with a real scan, so this shows where the bias comes from, not that a tissue-only dictionary would work on real data. Berman didn't test a tissue-only dictionary. Second, even 3B0's tissue signal isn't perfect: in Berman's Figure 7 its tissue error is about 14 percent, against about 3 for the good 2D recipe.

What this means for us depends on what our dictionary includes. Jack's notebook calls a tissue-only dictionary the manuscript's design, but the MrVox runs had the blood signal switched on, with ExtraVascOnly set to 0. If our dictionary is tissue-only, 3B0 matters much less; if it includes blood, it matters a lot. That's one of my questions at the end.

## Slide 20 · Scaling up: why the cluster

Everything so far ran on a laptop. A real dictionary is bigger, and that's where the cluster comes in.

One 3D simulation takes about 35 seconds on one core at 3 percent blood volume; more vessels take longer. A first real dictionary, say 25 radii, 15 blood volumes and 3 random voxels each, is about 1,100 runs. That's 10 to 12 hours on one core if every entry cost the same, and more once the high blood volumes are in. But each run is independent, so as 30 jobs side by side on the cluster it takes about an hour.

Two things make this manageable. Oxygenation is free, because one run can be rescaled to any SO2. And each run uses one core, so it maps onto array jobs the same way the MrVox jobs are already split.

My next step there is to get my account, test one job interactively, and then run a small array with Jack's template.

## Slide 21 · How I would turn this into my rotation

So here's how I'd turn this into my rotation.

The question is in two halves. First, does the way we simulate the dictionary change the blood volume, radius and oxygenation we report from GESFIDE? Second, once the dictionary is right, can a neural network replace the matching step? The second half is what's in my rotation agreement. The first half is what I'd add, because the network can only be as good as the dictionary it learns from.

Part one, over the next two weeks: build a GESFIDE dictionary the recommended way. Small on my laptop first, then full size on the cluster. It's done when it reproduces Jack's reference.

Part two, about three weeks: build a 3B0 dictionary with the same settings and measure the difference in all three parameters, with noise, and then on one lab dataset.

Part three, if there's time: a small network against matching.

The last week is the summary, the lab-meeting talk, and a clean repo. Even if part three doesn't happen, parts one and two are a complete result.

## Slide 22 · Where I would like your input

I'd really like your input on a few things.

First, is "does the simulation recipe change our estimates" a good first question, before the machine-learning part? Or would you rather I go straight to the network?

Second, what would help the lab most: a new dictionary you'd actually adopt, or a measured bias for the current one?

Third, which dataset should I use, and was it acquired with or without contrast agent? That changes what has to go into the dictionary.

And fourth, does our dictionary include the blood signal, or only the tissue signal? Jack's notebook calls a tissue-only dictionary the manuscript's design, but the MrVox runs had the blood signal switched on, with ExtraVascOnly set to 0. After this week's test, that's the question that decides how much the 3B0 problem matters for us.

## Slide 23 · Thank you

That's it. Thank you, and I'm happy to take questions.
