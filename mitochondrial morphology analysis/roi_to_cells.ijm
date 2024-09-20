// Prompt for output folder
output_folder = getDirectory("Select Output Folder");


// Get the title of the original image for later reference
originalTitle = getTitle();

length = 16 //length of selection


// Duplicate each ROI into a new image and save as TIFF
for (i = 0; i < length; i++) {
    // Switch focus back to the original image
    selectWindow(originalTitle);
    
    // Select the ROI from the ROI Manager
    roiManager("Select", i);
    
    // Duplicate the ROI
    //run("Duplicate...", " ");
    run("Duplicate...", "duplicate");
    
    // Get the title of the duplicated image
    duplicatedTitle = getTitle();
    
    // Save the duplicated image as TIFF
    saveAs("Tiff", output_folder + originalTitle + "_" + i + ".tif");
    
    // Close the duplicated image to free up memory
    close();
}

// Optional: Close the ROI Manager
// roiManager("Close");





