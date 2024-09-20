// Prompt the user to select a folder containing the z-stacks to be processed
input_folder = getDirectory("Select the Input Folder");

// Create a new folder to save the processed z-stacks
output_folder = input_folder + "Cropped_Stacks/";
File.makeDirectory(output_folder);

// Get a list of all files in the folder
list = getFileList(input_folder);

// Loop through each file in the folder
for (i = 0; i < list.length; i++) {
    // Check if the file is a TIFF image
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".tiff")) {
        // Open the current z-stack
        open(input_folder + list[i]);
        
        // Get the current z-stack title for saving later
        title = getTitle();
        
        // Make a centered rectangle selection
        // The original size is 273x293, and we are cropping 10 pixels from each edge
        makeRectangle(10, 10, 253, 273);
        
        // Crop the selection across all slices
        //run("Duplicate...", " ");
    	run("Duplicate...", "duplicate");
        
        // Save the cropped z-stack
        saveAs("Tiff", output_folder + title);
        
        // Close the current z-stack and the cropped one to free up memory
        close("*");
    }
}
