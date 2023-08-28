# RPE-Modeling | v0.0.2
Original code for models in "Internal timing-related dopaminergic dynamics are explained by reward-prediction errors" by Allison E Hamilos and John A Assad, biorxiv, 2023

> Credits
> 
>   Temporal Difference learning under state uncertainty with feedback:  Value function learning code was adapted from original code by John Mikhael (@jgmikhael). Dr. Mikhael's original code is included in the JMoriginalcode folder for comparison (<a href="https://github.com/jgmikhael/ramping">his latest GitHub version here</a>). Original derivation of this method is reported in <a href="https://www.cell.com/current-biology/pdf/S0960-9822(22)00036-7.pdf">Mikhael et al., 2022, Current Biology 32, 1077–1087</a>
>
>   Application of this method in models of the self-timed movement and perceptual timing tasks:  Allison E. Hamilos (@harvardschoolofmouse)
> 

## Dependencies

Runs on MATLAB 2023a+

## Getting started:

Running obj = CLASS_value_landscaping_obj will produce Figures 2 and 3 from the biorxiv manuscript. That's all you need to do! The uncertainty kernels will also display. 

### Simulated internal timing states and value/RPE landscapes
<img width="30%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/5477780d-c597-4191-aa1a-65142770361b"> 
<br>
<img width="40%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/8a073b6f-8b93-4497-94b4-1958e890013f">
<br>
<img width="40%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/9c02d581-09f3-45b2-8641-4f30eeda956f">


### Value function learning (the state uncertainty + feedback model is shown in the paper)
<img width="25%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/66dbc37a-6b0d-4684-9c39-1aa51d9c5b63">

<img width="25%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/9247e6cf-66ae-4bd5-b936-90a1e1c4a99a">

<img width="25%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/dd7660a0-328d-4316-ab6f-0234de591ae1">


<img width="25%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/45a28bf7-9b8d-483b-9b42-4a0a7c211b86">





