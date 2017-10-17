Andrew Berger  
October 17, 2017

This project contains Python code to process FMR spectra collected with a VNA.

To execute the program:
<ol type="1">
<li> Execute VNAFMR_GUI.py. This brings up the user interface window. </li>
<li> Enter file and processing information into interface window.
    <ol type="a">
    <li> You can enter the Filename manually, or using the "Select File" button. 
       The full-path must be included </li>
    <li> Files are expected to be tab-delimited, with no column headings </li>
    <li> External field units are expected in Oe </li>  
    <li> Filenames should be of the form:  
        "directory/file description X.XXX GHz Final etc"  
        where:  
        <ul>
            <li> the frequency is included (e.g. X.XXX = 10.000) </li>
            <li> " GHz" demarcates the end of the frequency digits (a space is 
            necessary between the digits and the "GHz" label) </li>
            <li> "Final" is included somewhere in the filename to indicate that the given 
            file includes the final spectrum (after repeated averages, e.g.) </li>
        </ul>
        </li>
    <li> "Column Selection" indicates the data columns for field-swept VNA-FMR spectrum
        <ul>
            <li> Field </li>
            <li> Re(S_21) </li>
            <li> Im(S_21)  </li>
        </ul>
        </li>
    <li> "Frequency Selection" sets the frequency ranges for fitting.
        <ul>
        <li> "Spectral Fits" sets the range for the field-swept spectra to be fit. 
        If the "Minimum Frequency" and "Maximum Frequency"
        are set above and below the minimum and maximum frequencies that were measured, 
        all spectra will be fitted. </li>
        <li> "Linear Fits" sets the range for Kittel, Damping, and Inductance 
        linear fits vs. frequency. This can be a subset of the results from the 
        spectral fitting, in case the resonance field, linewidth, or inductance
        vs. frequency becomes non-linear. </li>
        </ul>
        </li>
    <li> "Automatic Windowing" sets the field range over which the spectra are fitted.
        1.5 linewidths has provided good fit results with our experiments. Linewidths
        are given in units of the linewidth defined in the Polder Susceptibility
        (see "PythonVNA_FMR.pdf"). The full fitted range is given by the fitted 
        center resonance field +/- (# of linewidths)*dH. Therefore, for a user-input 
        value of Window Size = 1.5, the total fitted field range will be 3*dH. </li>
    <li> Results are reported in the same directory as the source data, in a folder 
        specified by the user in "Output Folder Name". </li>
    </ol>
    </li>
<li> When "Run Fits" is clicked, the Python console or shell executing the code will 
    print a progress report including the currently fitted frequency, and the field 
    range to be fitted (in Oe and in # of linewidths). </li>
<li> When fits are completed, a report will print to the Python console or shell
    including the results of the Kittel fit, Damping fit, and inductance fits 
    vs. frequency. </li>
</ol>
