This README Accompanies:
------------------------

*Human scent guides mosquito thermotaxis and host selection under naturalistic conditions*

Diego Giraldo*, Stephanie Rankin-Turner*, Abel Corver, Genevieve M. Tauxe, Anne L. Gao, Dorian M. Jackson, Limonty Simubali, Christopher Book, Jennifer C. Stevenson, Philip E. Thuma, Rajiv C. McCoy, Andrew Gordus, Monicah M. Mburu, Edgar Simulundu and Conor J. McMeniman

--------------------
Corresponding Author:
---------------------

Conor McMeniman
cmcmeni1@jhu.edu

------------------------------
Files in this directory:
------------------------------

- MosquitoFrameReader.py: Process raw recording files and perform background subtraction and ROI detection.
- MosquitoTrajectories.py: Form trajectories/detect landings by linking mosquito occurrences across frames, using various heuristic rules to filter out non-landing events.
- MosquitoTrajectoriesExport.py: Export detected landing trajectories to CSV for subsequent analysis.
