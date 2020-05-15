// Arpag, Lawrence et al., PNAS, 2020
// Macro for making multiple kymographs from a two-color time-series.
// This macro first makes a max projection (chose "No" for "Use old max projection?" question)
// The initial time series needs to be split into channels while running the macro for the first time (chose "Yes" for "Split Channels?" question)
// User then needs to draw segmented lines on microtubules and add to ROI manager.
// Once the lines were drawn on desired microtubules, user needs to click "OK" to proceed.
// Kymographs are made using built-in "Multi Kymograph", then combined into a stack.
// If user wants to draw more lines, the macro can be rerun while ROI manager keeps the old ones (or previously saved ROIs can be loaded).
// Before rerunning the macro, user needs to click on the initial time-series window.
// If the max projection window is not closed, user can chose "Yes" to "Use old max projection?" question.
// If split channel windows were not closed, user can chose "No" to "Split channels?" question.

initialname=getTitle;
selectWindow(initialname);

choices = newArray("Yes", "No");
var choice ="";
Dialog.create("Select choice");
Dialog.addMessage("Use old max projection?");
Dialog.addChoice("-->", choices);
Dialog.show();
choice = Dialog.getChoice();
if  (choice == "Yes") {
	print("using old max projection");
	selectWindow("MAXProj-"+initialname);
}
else{
	print("making new max projection");
	selectWindow(initialname);
	run("Z Project...", "projection=[Max Intensity]");
	rename("MAXProj-"+initialname);
}

choices = newArray("No", "Yes");
var choice ="";
Dialog.create("Select choice");
Dialog.addMessage("Split Channels?");
Dialog.addChoice("-->", choices);
Dialog.show();
choice = Dialog.getChoice();
if  (choice == "No") {
	print("using already splitted channels");
	selectWindow("MAXProj-"+initialname);
}
else{
	print("Splitting channels");
	selectWindow(initialname);
	run("Duplicate...", "title=dup duplicate");
	selectWindow("dup");
	run("Split Channels");
	selectWindow("C1-dup");
	rename("C1-"+initialname);
	selectWindow("C2-dup");
	rename("C2-"+initialname);
	selectWindow("MAXProj-"+initialname);
}



setTool("polyline");
run("Select None");

roiManager("Show All");
roiManager("Show All with labels");


beep();
waitForUser("action required", "Draw segmented lines on microtubules.\nAdd each line to ROI manager.\nClick OK to make kymographs.");

nROIs = roiManager("count");
print(nROIs);

// make kymographs

for ( i=0 ; i < nROIs; i++) {
	// make kymograph with first channel
	selectWindow("C1-"+initialname);
	roiManager("Select", i);
	run("Multi Kymograph", "linewidth=3");
	rename("C1-Kymograph-"+i);

	// make kymograph with second channel
	selectWindow("C2-"+initialname);
	roiManager("Select", i);
	run("Multi Kymograph", "linewidth=3");
	rename("C2-Kymograph-"+i);

}

run("Images to Stack", "method=[Copy (top-left)] name=Stack title=C1-Kymograph-");
rename("Stack-C1");
run("Images to Stack", "method=[Copy (top-left)] name=Stack title=C2-Kymograph-");
rename("Stack-C2");
// chose appropriate channel colors based on the raw data. c1=red, c2=green, c6=magenta etc.
run("Merge Channels...", "c2=Stack-C1 c6=Stack-C2 create");

rename("Stack-"+initialname);
