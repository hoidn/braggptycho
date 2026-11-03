<meeting 1, anl>
direct beam, saxs profile
beam path direction
if standing at saxs beam line, 
in plane peak: projection back of the crystal
bragg peak above the plane: vertical (horizontal) segment

if translate beam, front and back surfaces don’t move (ie can’t tolerate thickness variations). has the benefit of mathematical simplicity. 
(nothing about bragging geometry, solve objections 

forward problem: ray path is invariant, solve for direct beam view
direct beam view: all details of diffractometer and sample

uniform 2d slabs, 3d structure (3rd dimension will matter. we want to use the 2d saxs algorithms… fit the information in a 2d array that makes the least error with respect to your diffraction. walking in this direction one wants to pick samples… membrane of silicon with some dislocations, essentially a 2d problem). (things are observed at some angle, the pixels don’t know about this)
another paper: with single angle, solve for variation in depth direction (but it’s really fragile) 

30 pixels from the center: down to 10s of photons per pixel. 
another difference vs saxs ptychography is that there’s a rocking curve for each bragg peak.
mostly fix at one bragg angle and collect additional information. question: where on the rocking curve do you want to be?

040 phi angle of 45 deg incident to beam. 45 deg inplane det displacement. the outcome: kb focused beam (should be a square; at a 45 deg with beam transiting through the sample it’s a rectangle with a band through the middle). 
if you figure out a small ray path (as if looking at beam from 45 45 pov) 
do a convolution of two arrays to get the effective beam. this is how to distill the geometry
phi0: beam normal to sample. one motion for sample rotation. two in plane diffraction peaks, one mixed in / out of plane. convolution kernel changes as you go. the elongation falls out 
dislocations induce a deflection. in the dislocation field it gets interesting bc in 3rd scan no dislocation is visible (zero dot product)
unified problems where the target isn’t phase and amplitude but the u field that reconciles all bragg peaks and positions. (unknown burgers number and ..)  different aspects of the crystal
depth resolved or depth averaged?: a mix, bc it depends; you have some amount of tomography; 3 different projections; observation angle moves with the detector, so you have to triangulate. 
at high energy beam lines (access many peaks) this triangulation works well. but typically you don’t have the energy to reach high q / multiple scan all the raster peaks. here: trying to save time by taking just dataset on the bragg peak. (split peak -> need multiple phi angles; works against you if you need multiple bragg projections). 
one model connects them; how do you get the lattice displacement field that connects them.
steve: would need detector registration. 
david: defined titanium e ben deposited; different shapes.
second image as a marker that you can use to align multiple images. fluorescence marker, might create problems but allows you to align the images per pixel. 

membrane wants to strain / stress. Maybe there’s feature recognition stuff that’s very automatic. 
1st goal: reconstruction of individual bragg peak, then figure out ho to combine them and deal with combining when there are offsets. if you find that this isn’t working
if field of view has curvature bragg peak moves inside rectangular aperture, then everything is moving bc of the lattice strain; then …. how much in phi are you off in the bragg peak bc to explain a bragg peak down here in terms of phase you start to have to pile up phase ramps, and this becomes difficult to interpret. 
if working with domains with stress at grain boundaries start to have problems; there are things you can do; this goes into the 3d model; could be that this 2d scan works, but we’re definitely not there. 

when is amplitude / phase relevant vs. displacement field. to go directly to displacement you’d have to solve for all of the offsets. this has been worked out in simulation (?what’s the reference?) 
incident beam -> dislocation -> off of surfaces on your own. assume no dynamical diffraction. having not enough photons helps. if you have strong bragg peak you’ll have … extinction factors, these issues will be ignored for now. 

how do you deal with sample geom non perpendicular to beam? do you see one edge of the sample at different focal.? : at 26. you’ll see it but at max 40 you won’t bc you have sample scan instead of beam scan; beam hitting the surface at always the same relative point; but in general you’re right; if you’re at an oblique angle and move the beam you’ll be off focus. …. numerical aperture, shortest steps of focus.
why? resolve features smaller than the beam. 5 nm beam. what are you scattering from that bragg ptycho will solve for you? unclear
would prefer to use modes focusing ? …. bit circular bc what you want … longer focal length lenses, ability to itnerfere over many periods. 
wavefront is curved wave; 

what are the issues? sensitivity of the inverse problem, initialization? if you work with materials with dislocation stresses etc., 2d models don’t work. here, the origin moves with the crystal planes. (any bending means that you’re chasing this origin). 
becomes very hard to figure out how to do this without multiple phi slices. bragg peak is somewhere but most of reciprocal space is empty, cut you still have to take the time to scan to see where the peak went. 
in simulation: start with slab geometry and then start introducing defects.        

get the 111 data from swan. use the recon that’s working and see if we can get it. 

</meeting 1, anl> 

<meeting 2, slac>
Speaker A: Um, is that— or like a partial dislocation. So calling sometimes to see if you can visualize images of the dislocation. So one thing that's interesting is that—

Speaker B: No, 20 is not the only gap.

Speaker A: It's 20, you're like— 5, 8, and 9.

Speaker B: It's only 9.

Speaker A: You're probably right.

Speaker B: Like problem layer by layer.

Speaker C: Yeah, so you can say that's what we think. So if you look at this type of photo, that it might be twin columns that are twisting over, look at this thing. We don't think, but that's something we need to see.

Speaker A: But one thing that's interesting, so like these 2 peaks are from, let's say, twin A, and then this peak is from twin B. And you can see that in twin A, the mosaic spread is in this direction, and for B, the mosaic spread is in this direction. And they don't even line up. So in the 6-nanometer case, it's not a— it's like in one machine, it's like this. In one direction, it's like 176, and the other direction is 184 degrees. So there's some sort of Symmetry breaking. It's not like we have one domain that you can split in both directions. Like, one domain is only going to the left, the other one is only going to the right. And then there's also a monoclinic distortion that I need to quantify, but it exists in the 6-nanometer film, and that distortion actually increases as we go to the 40-nanometer film, but we lose this— they now all are dimensionally So the domains become 180 degrees separated past this angular formation, but then we get all of this crazy diffuse scattering. It's like kaleidoscope looking. We see there are 3 main diffuse scattering that we see. We see this like phonon-looking scattering. We see these rods, which I think are like bunches of stacking faults. But actually quantify that. And then we see also partial work. So we have basically a situation where something gets strained in the domain, and then that breaks after a certain point, causing all of these defects. And we think this is— we can find the mechanism of the polymer and switch to polymer by understanding other polymers. So we took the table to At Argonne to try and figure out—

Speaker C: You're on top.

Speaker A: Oh, uh, yeah, thanks. Um, so we had, uh, 2 zone plates that we used. This is a little weird because they sort of just didn't have the right zone plate.

Speaker C: They were having issues.

Speaker A: We used 2 different zone plates. One was like really weak density, but it had a much higher focus, so it was like point managers. And the other one was about— I need to get the specs for installed, but it's like 150 nanometers plus. And it was like way brighter, so we could get a ton of speckle.

Speaker B: Where were you operating? In both.

Speaker A: We were operating in both. Accessible data that's kind of like HyperTree ready if we want to do that. Yeah, we took like a series of focus positions. So we have, yeah, we have In 2020.

Speaker B: Oliver, I see you're here now.

Speaker C: Can you hear us? Yes, we are. Oliver, can you hear us?

Speaker E: Yeah, I can hear you. You went out for a little bit, but you're back.

Speaker C: Okay, um, maybe we should back up for a little. You are just showing some data, but this whole thing— can you go to the previous slide? You're looking at first, uh, sort of make it very clear, we are looking at a 6-nanometer, uh, coordinate, uh, a full data, SSR data enhanced structure that we were talking about, uh, and We collected data on ID turning with a sun plate and maybe 2 sun plates and several different settings of the sun plate. We are on in-focus. It's something sometimes we are, you know, the sample is behind the focus. And David was just going to show some data that we collected on You can say that.

Speaker A: Yeah, we have 6 nanometers.

Speaker C: We have basically—

Speaker A: there's like a transition in growth mode as you cross 10 nanometers. So we have a 6 nanometer film where the growth mode is still layer by layer. We have a— not shown here, but we have an 11-ish nanometer film that's just past this switch. And then we have a 40 nanometer film that's well past the switch and has basically a creation of these huge column or grains. We also have a gold sample, and that's interesting for like a more methods standard sample.

Speaker C: Do you have any questions, Oliver, or should we go into data?

Speaker E: Um, no, you can go into data.

Speaker A: And actually on this slide, let me say one more thing. So we took sort of 3 different darkfield maps, one of the specular diffraction, which is the one-on-one. And then we found, for example, 2 off-specular peaks, one for each of the domains. So one was the 202 for one of the domains and then the 131 for the other domain. So they're like the 2 twins. So we basically find both of them to make the maps.

Speaker C: The spectra, I don't know if you could see that they were twins.

Speaker A: It was like Remarkably flat, the intensity. There was still a lot of speckle there, but in the intensity there was almost nothing.

Speaker C: But in the all speckle, yeah.

Speaker A: So this is a gold—

Speaker C: okay, so it's simple.

Speaker A: So we have a gold sort of standard film. This film is like a really super pristine MBE film. But it still has this twin domain structure.

Speaker C: Can you see online?

Speaker E: Yep, I can see.

Speaker A: Okay, so we have this sample, and this sample is like acting as a standard because it has the same kind of domain structure, but The piece of order, and then, and I mean, like, you can, you can sort of like fit all of this to theory and kind of make it seem as though it's— The second ring, because it's with its own center, it has a shadow. It doesn't completely get the third order, uh, completely removed. It's not perfectly Um, so for this sample is sample 1, but we could— we had like, these are like kind of what the domains look like for an offsite. This is actually just some text here.

Speaker C: Absolutely.

Speaker A: But like, just so like the cold zone is like a good standard because, you know, you can see we have like, you know, a million versus in this, uh, much larger domain.

Speaker E: And so forth.

Speaker A: But this is not necessarily the— that's like the science of it is not.

Speaker C: But we should try to maybe think we can reconstruct the goal as following by the spike potentials and by opinion samples. Yeah. And then this is the place where I was telling you guys that Can we use intensity maps or some of the correlation methods we were thinking of? There are 2 kinds of correlation methods we were thinking about. Can we use them to create— what does this tell us? This wave structure, is this related to some image? Can we use that to encode some starting space for solving the reconstruction? Yeah. You want to show a couple of diffraction patterns and then we can talk about it?

Speaker E: Yeah.

Speaker C: So this is—

Speaker A: Yeah, well, let's say, so this is— let me start with this. So this is the finished sample.

Speaker C: This is what in general what the specular peaks look like.

Speaker A: So there's quite a bit of speckle, and they have— they're all— they're typically all elongated like this. We can find that as we raster across the sample, these speckles are Physically moving across the detectors, so we can see that they're actually being tracked as we move the stand. As we raster the probe by like 100 nanometers, we don't know that we can say it looks like they're moving along. Well, I should say, if you look at one spot, it kind of continuously moving and—

Speaker C: MR spot movement. This is, uh, this one is probably—

Speaker A: this image is focused, but we can view this.

Speaker C: So what, what does power affect it?

Speaker A: Well, so all of these were taken with, uh, like 100 feet with a 20-nanometer spot size mode plate. You didn't see much speckle. It was just, uh, all of the speckle, all of the patterns were Psychometrics is a 150 nanometer spot.

Speaker C: So one initial thought we had is even if you can't just look at this, cycles have this orientation, and if you can look at 2 diffraction patterns next to each other, then you find correlations. You can actually see firstly just a raw correlation of the other 2 images. Pixel by pixel, and if they're moving, we can have another— can do something even more smarter by saying, you know, sort of rotate this thing so the lines also parallel and perpendicular to that. And can we actually come from evidence saying that, you know, vertically along this, the correlation is changing on average at this rate, and on this other, it may be changing at that rate, and come up with a line fit. And then we do this for the left image and then come to the right and then move to the top. And so we basically have 4 numbers for each image and then we create those 4 numbers and drag them across the thing and see if that gives us some sort of a substitution. Probably not a real, but if we have something like this, only what I was thinking, that might Then we can sort of see if that— how it connects to the actual sculpture. We can then connect that to the center of mass or intensity maps. And if they are similar, then maybe we can use this to see the reconstruction. Yeah.

Speaker D: Right.

Speaker C: Another thing we saw also in the specular, we actually saw satellite peaks here. So these have some—

Speaker A: Yeah, so you can see it very clearly here, but if I like sum over an entire map, then you can see these ones as well. much clearer. But there's some like background speckle, and then we also— You can look at maps, right?

Speaker C: Is that how to find those satellite features?

Speaker A: Yeah, so here's an example here.

Speaker C: Here's a map.

Speaker A: Oh yeah, here it is. So this is like how I— you pick out one satellite feature and you can see that map here. So this is telling us that there's some periodicity that's like, uh, 25 nanometers or so, um, that we're just not able to pick up with those probes because this probe is 100 nanometers. So we don't have the resolution to actually just see it to be intense enough. Um, but this is, you know, something pretty cool. Pretty well-ordered period there at 25 nanometers or so. And then as we grow thicker, so this is a 12-nanometer sample, this is just past the critical thickness. We don't see—

Speaker C: So say what that is.

Speaker A: This is an intensity map of another offset region.

Speaker C: But you can see that it's nice order.

Speaker A: So we don't see these satellite peaks, but now we see a much larger period. So these are like— I think the step size of each pixel is probably like 2 nanometers, so this is a couple hundred nanometers period. And actually, as we zoom out, we see there's like a second order in the other direction. This is what the speckle looks like. And what's happening is that— What's happening is that like those streaks going in this direction is because These patterns actually have like a break in it. So you can sort of see that there's— zoomed in, there's another pattern going in this direction that then becomes much more apparent when we zoom out.

Speaker D: Do you have a general, like, how far in angle space those 2 are separated? Like the angle between those 2 directions?

Speaker A: They're not 90 degrees, but they're a little bit above 90 degrees. I very roughly by eye just had to pull that out, but I— yeah, it's not nice. So yeah, so basically once you grow past 6 nanometers, you don't have this really tight periodicity. We have a much larger grain scale that we can actually see with a nanoscope. And then when we go to the thickest one, We still see kind of a stripe pattern. There's something here, uh, it's— but it's like much harder, it's much weaker. Um, that was taken with a smaller zone plate, then we moved to the larger one so that we— we sort of realized that we should move to the larger one because we're getting a lot more speckle information. The speckles are kind of a lot smaller. extent in the detector. With the large zone plates, this is the sort of maps that we're taking. So we still see now the stripe pattern has like way larger periodicity, and there's like— they're much more wiggly, they're not as well ordered. So this is kind of what the data looks like, but in general we have the same sorts of maps for every single one where we have off-spectacular reflections for 2 domains. We also have like rocking through the actual peak, so you have like several angles through the peak. So we have kind of like an RSM of the speckle at every point, not like a densely measured RSM, but, you know, 3 angles through the peak. So this is, this is the kind of data.

Speaker C: Question. On your repeat.

Speaker E: So what would be the most useful dataset for us to start with if we're looking at starting to develop?

Speaker C: Someplace would be covered, I think, because we know what the tech domains look like. The most interesting is this.

Speaker B: Yeah, I think the gold does not have the photo.

Speaker A: The from the light in the photo.

Speaker C: Oh, it did have the yeah, it did have a shadow from light. Make it three make it. So we kind of know the domain structure of the gold looks like. So the gold yeah, do we know?

Speaker A: We have like EBSD of the domains from sort of sister samples. We know like what length scale, what kind of shapes the domains are in general. So that's why we have it as a kind of analog.

Speaker C: They're not structured like this.

Speaker A: They're not structured. They're wiggly and different shapes and kind of random.

Speaker C: So in that sense, looking at symmetry is always a good question. Good, or maybe the example are the ones to begin with, right?

Speaker A: Yeah.

Speaker C: So there are cycles that are defined. So maybe the thing to play around with, or however you want to do this, go to the steps of the example. Maybe that's the thing. Yeah, if you look at data like this, maybe not even do reconstruction, but what are our thoughts on looking at, you know, image by image? Is there some image correlation? What do we get from that? It's not really dichotomous equal mitigating, it's dichotomous information. Yeah, step sizes were 400 nanometers.

Speaker A: We have a bunch of datasets ranging from like 20 nanometers to 100 nanometers step size.

Speaker C: You can do a correlation map and that works. Then we get some sort of an image, if you figure out how to plot it. With a 10-nanometer resolution, and that might be a good starting point. It may not tell us exactly what it is, but it might give us some idea of what the features look like. So maybe that's one place to start.

Speaker B: I have decent probe information from the beamline.

Speaker A: We can get it.

Speaker C: They do have it. All right, so we can get a probe also, but Maybe the simplest thing to start is maybe just this, right? Can you look at image by image correlation? Maybe all of you have some thoughts. If I have this bunch of images, how do I find correlations across? How do I plot positional correlations from one point to the next point? I was also thinking taking, you know, one point and looking at image right off it and on top of it or whatever and maybe getting some correlations. Maybe initially just measuring, you know, uh, a contrast between 2 images, but maybe we can do something more complicated, uh, if there are speckles that move in different directions or some tracking or some similar— I don't know, any thoughts that you have on that? Maybe that'd be the kind of cool thing if we can figure out what we can do with the correlations. that might actually open up scientifically a lot of things we can do. We cannot do reconstruction. We can do this and, you know, reason out that we would not otherwise get. But I don't know, what do you think?

Speaker E: Yeah, I mean, if we start with simple things like cross-correlation There might be a second Fisher information approach.

Speaker A: What is Fisher information? What is Fisher? I don't know what that is.

Speaker E: It's like an information— oh, go ahead. Oh, I was going to say it's a measure of mutual information from more of a Bayesian or information theoretic perspective. point of view.

Speaker C: Okay, okay, okay.

Speaker A: Um, so, so that would give us like a scalar metric kind of.

Speaker C: Um, how can L all the way go back to the base images? So I'm on a hard drive, but we also This endpoint here, we put it on global support here. Can we put it on S3DF or something?

Speaker A: Yeah, we— on MPC, we don't have a space.

Speaker C: However, would we want to have the data?

Speaker E: Yeah, it'd be nice to play around with this, definitely. Um, S3DF is convenient, but Anything that works, really.

Speaker C: You have a location on SPCS?

Speaker E: How big is the dataset?

Speaker C: It's actually—

Speaker A: each map is pretty heavily compressed, so it's only a few, like 5 to 10 gigabytes, I think, per map. And then uncompressed is like a few hundred gigabytes. But so the data itself, we have, I think, 300 gigabytes of data uncompressed.

Speaker C: So each file is one map, and one map is— Yeah, maybe 100 gigabytes.

Speaker A: Uncompressed is, yeah, a few hundred gigabytes.

Speaker C: What do you think?

Speaker E: Um, could— does Google Drive work? Do we have enough, um, allocation on there?

Speaker A: Why don't you, or even like a OneDrive link? Are you thinking for just having the data?

Speaker E: Yeah, just for sure.

Speaker C: Why don't you, why don't you just talk to all the— Yeah, we can do that.

Speaker A: I can also just give you a Globus link, like for the Argonne, the whole Globus endpoint link, and then if that's convenient, or whatever's convenient for you.

Speaker E: Oh yeah, that's actually convenient, like in Globus, yeah.

Speaker C: Okay, okay. And then maybe, um, maybe we can point you to the datasets to look at because there's probably a lot of labels. Maybe pick one that— maybe this one is a good one.

Speaker A: Yeah, we have a lot because we weren't exactly sure. Like, we could see that the speckles were moving, we could kind of see over what length scales we We can still track them, but we weren't sure like what is the best data. So we have like a lot of datasets with different step sizes and different defocuses, and so we have a lot. So we can chat about what's the best one to start at.

Speaker D: Yeah, make sure you have the link visible.

Speaker A: Yes, yes, I'll share with you this, and then we have sort of every single scan listed here. So like basically all the details about motor positions and everything, the intensity maps. So I can share with you.

Speaker E: That sounds good.

Speaker C: We should meet regarding anything? Every week? Every 2 weeks?

Speaker E: I'm sorry, every week or every other week?

Speaker A: Either is fine for me. I think I would be fine with it. I would be fine with this one because Too frequent?

Speaker E: I think I can do every other.

Speaker D: Does Oliver have time to do something?

Speaker C: Oliver, what about you, bud?

Speaker E: I think every other would work for now. I'll probably have more bandwidth starting in the fall.

Speaker C: Yeah, yeah, that works.

Speaker E: Like September.

Speaker A: Cool, cool.

Speaker C: And Friday around this time work?

Speaker E: Yeah, that works. Yeah.

Speaker C: All of us should be involved.

Speaker E: Oh, sorry, my audio is breaking in and out.

Speaker C: I didn't I was saying, should we involve Albert and Steve, or do we just keep it to Slack right now?

Speaker E: I think Steve especially would be interested in this kind of thing, is my impression.

Speaker C: Yeah, I think so too.

Speaker A: So when I set up a meeting, I'll talk to Steve.

Speaker C: He'll ask me to join. Cool.

Speaker E: Cool. Okay, sounds good.

Speaker C: मेंशन तो है इस बार, एंड इंटरचेंजेस क्योंकि वो सेंटर से आईडी पॉइंटेड, तो दूसरा ही है कि मतलब सम, ये है तो वो जो डेटा पाइपलाइन होती है वो है, यार, आउट ऑफ़ मीटिंग विथ स्टीव एस वेल, तो डॉक्टर डॉक्टर, ओके, मैं यू हैव एनी क्वेश्चन?

Speaker B: Nothing too specific. I mean, there's going to be some points. The correlation thing is critical.

Speaker C: Yeah, it was. It was supposed to be done.

Speaker B: I mean, I don't see much work.

Speaker C: We can find something that compares. We'll find something. Yeah, but since it's powered by a book, the book will definitely give us something useful. It will throw us some information. What do they do? Probably a hard question I'm going to ask you, Matt, but how is it just a pattern pattern for a researcher?

Speaker B: I mean, I think it's actually fairly the same, more that it's complicated.

Speaker C: Okay. So yeah, lots of mini-masters you can follow and do. Right. So this collapse is not very—

Speaker B: it's not going to result in a simple—

Speaker C: that's what I'm wondering.

Speaker B: No, I know it's recognizable.

Speaker C: So I'm wondering if you can get a collision map or some sort of map that says this is approximately what the real space would look like. Yeah, it's just—

Speaker B: I mean, that's going to be by definition limited to the resolution of Test size or actually the size of the probe. So then to try to look at finer features, you need a small probe.

Speaker C: So there's a spectrum between low-level structure and high-level structure. Fine subpages and high-level structure. Yeah, that'd be an interesting thing for Steve and you to think about a little bit.

Speaker B: I guess the scanning is probably always along the sample surface.

Speaker A: Yes, so I guess we—

Speaker B: like, you're scanning parallel to the sample surface?

Speaker C: No, no, like I said at the start, I don't have the scanning—

Speaker B: Oh, you're scanning the zone plane then?

Speaker C: I'm in the zone plane, so at the start—

Speaker B: Okay, so actually the— Distance from some place to sample is changing from one step to another.

Speaker C: That's actually kind of complicated. It's a lot.

Speaker A: I guess for a small field of view, it's like 15 degrees into the navel.

Speaker C: Good place to take this.

Speaker A: But it's only 15 degrees.

Speaker C: Yeah, I thought I might not want to do it. It should be like a constant flow too, once you figure out what that flow is. No, but the, the flow first.

Speaker B: Yeah, the flow.

Speaker C: You know, the flow. So that's—

Speaker A: yeah, the whole— quite a— yeah, fluency. Yeah, I guess hypothetically, like, this is pretty—

Speaker B: I don't know, we could have had like last week.

Speaker A: Like, you could actually.

Speaker B: Yeah, yeah.

Speaker C: So that would—

Speaker A: these permissions, the magic way is if you can scan a sample.

Speaker C: I understand that way. That's why I think the thing is one place that's pretty worth doing it. Yeah, yeah.

Speaker B: Like, if you could scan a dump list along those Yeah.

Speaker C: All right. Wait, let's wait in two weeks. I'll send out an email to Steve, Ryan, Joe, put all of us on this.

Speaker E: Yeah.

Speaker C: Put out your agree or whatever works. Yeah. Thank you. Two weeks.

Speaker D: Okay.

Speaker B: I will be.

Speaker A: I'll have to make the meeting. Online, that's fine.

Speaker C: No, Oliver's online. Okay, great.

Speaker B: Thank you, Oliver.

Speaker A: Yeah, thank you guys.

Speaker E: Talk to you later.

Speaker C: Bye.

Speaker B: Bye.

</meeting 2, slac>

<followup email>
# nanoprobe data sharing from ID26 on NbN thin films

- **From:** "Garcia, David Andres" <davidagw@slac.stanford.edu>
- **To:** "Hoidn, Oliver" <ohoidn@slac.stanford.edu>
- **Date:** Mon, 17 Aug 2026 07:35:56 +0000

---

Hi Oliver,

I apologize for the massively delayed email about our nanoprobe data! If you have an Argonne domain account, you can access the data through the beamline globus endpoint using the following instructions: https://git.aps.anl.gov/DM/dm-docs/-/wikis/DM/HowTos/Getting-Data-From-Globus, with collection APS:DM:26ID and use the path /2026R2/S26_20260706/data/ as the source endpoint to start transfer.

If you cannot access that, let's figure out how I can transfer them to S3DF for you? I've attached a few files through OneDrive that are representative and probably good places to start looking at the data.

Here are the two notebooks from the experiment:

<link omitted; see links.local.md>

<link omitted; see links.local.md>

The scans 308 and 308 are of the (202) reflection of the 10 nm NbN thin film taken with 40 and 20 nm step sizes. The scans 312 and 313 are the same sort of scans but of the (131) reflection. The (202) and (131) are each from the different 180 degree domains. Each of these scans are taken at different (but nearby) regions of the sample because we had to limit beam damage. Each of these are taken with the zone plate in focus, which means the spot size is probably about 150 nm. I still need to get the probe from the beamline.

The hdf5 data structure is /entry/data/data, with each step as a 2d detector image chunked in data. The scans are rastered in a grid, with the steps defined by positioners ZPX and ZPY in /entry/scan/Positioner where ZPX is a grid of points scanned at each ZPY. The nth detector image in data corresponds with the nth value in ZPX.

- [scan_308_100308 2.h5](<link omitted; see links.local.md>)
- [scan_309_100309 2.h5](<link omitted; see links.local.md>)
- [scan_312_100312 2.h5](<link omitted; see links.local.md>)
- [scan_313_100313 2.h5](<link omitted; see links.local.md>)

Please let me know if you want to chat at all before friday about the data or have any questions,

David

</followup email>


<references>
../PtychoPINN
</references>
