filepath1 = File.openDialog("Select a File");
open(filepath1);
dir = File.getParent(filepath1);
fname = File.nameWithoutExtension();
fname2 = replace(fname,"_MMStack_Pos0.ome","")
waitForUser("select square qround neuron");
setTool("rectangle");
run("Align slices in stack...");
run("Enlarge...", "enlarge=60");
run("Crop");
saveAs("Tiff", dir + "/" + fname2 + "_reg");
//close();
waitForUser("select ROI around neuron");
run("Set Measurements...", "area mean modal min median stack display invert redirect=None decimal=3");
saveSettings;
setOption("Stack position", true);
       for (n=1; n<=nSlices; n++) {
          setSlice(n);
          run("Measure");
      }
restoreSettings;
saveAs("Results", dir + "/" + fname2 + "_neuron_results.csv");
run("Clear Results");
waitForUser("select background ROI");
saveSettings;
setOption("Stack position", true);
       for (n=1; n<=nSlices; n++) {
          setSlice(n);
          run("Measure");
      }
restoreSettings;
saveAs("Results", dir + "/" + fname2 + "_background_results.csv");
run("Clear Results");
close();