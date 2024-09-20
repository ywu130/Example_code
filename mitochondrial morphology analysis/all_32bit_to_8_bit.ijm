// Define the folder where the TIF files are located
//inputFolder = "/Users/yuanyuan/Downloads/200nM deepred/d10_n/d10_n_9-81/d10_n_9-81_crop_scale_L1/selected slices/dedonvolutionized/";
//outputFolder = "/Users/yuanyuan/Downloads/200nM deepred/d10_n/d10_n_9-81/d10_n_9-81_crop_scale_L1/selected slices/dedonvolutionized/8-bit-for-mito-analyzer/";

// Prompt the user to select the input folder
inputFolder = getDirectory("Select Input Folder");

// Prompt the user to select the output folder
//outputFolder = getDirectory("Select Output Folder");




// Get list of TIF files in the folder
list = getFileList(inputFolder);

// Create a new folder to save the processed z-stacks
outputFolder = inputFolder + "8_bit/";
File.makeDirectory(outputFolder);





for (i = 0; i < list.length; i++) {
    // Check if the file is a TIF file
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".TIF")|| endsWith(list[i], ".tiff")) {
        // Open the TIF file
        open(inputFolder + list[i]);
        
        // Get the title of the currently open image
        title = getTitle();
        
        // Convert to 8-bit
        run("8-bit");
        
        // Save the 8-bit TIF file
        saveAs("Tiff", outputFolder + title);
        
        // Close the currently open image
        close();
    }
}
