.onAttach <- function(lib, pkg){

  packageStartupMessage(
      "\n",
      "######################################################################################\n",
      "\n",
      "Have fun with GLAD\n",
      "\n",  
      "For smoothing it is possible to use either\n",
      "the AWS algorithm (Polzehl and Spokoiny, 2002,\n",
      "or the HaarSeg algorithm (Ben-Yaacov and Eldar, Bioinformatics,  2008,\n",
      "\n",
      "If you use the package with AWS, please cite:\n",
      "Hupe et al. (Bioinformatics, 2004, and Polzehl and Spokoiny (2002,\n",
      "\n",    
      "If you use the package with HaarSeg, please cite:\n",
      "Hupe et al. (Bioinformatics, 2004, and (Ben-Yaacov and Eldar, Bioinformatics, 2008,\n",
      "\n",
      "For fast computation it is recommanded to use\n",
      "the daglad function with smoothfunc=haarseg\n",
      "\n",
      "######################################################################################\n",    
      "\n",
      "New options are available in daglad: see help for details.",
      "\n")

}
