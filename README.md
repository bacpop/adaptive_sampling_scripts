# adaptive_sampling_scripts
Scripts for adaptive sampling analysis

analyse_RU.py is an edited version of the script available [here](https://github.com/SR-Martin/Adaptive-Sequencing-Analysis-Scripts).
This edited script allows removal of multimapping reads from [Minimap2](https://github.com/lh3/minimap2) and filtering based on % identity to the aligning reference.

split_by_channel.py can be used to split reads into two separated files depending on specified channel. For use when adaptive sampling has been used on some but not all channels.
