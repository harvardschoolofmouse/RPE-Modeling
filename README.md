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

<img width="40%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/8a073b6f-8b93-4497-94b4-1958e890013f">

<img width="40%" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/9c02d581-09f3-45b2-8641-4f30eeda956f">


### Value function learning (the state uncertainty + feedback model is shown in the paper)

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/4f8c63f5-c78f-4063-8eca-5520f8a851ee">

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/d20dc79a-a6d5-4d14-ab7c-8aa117e5d1b0">

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/78b1aaf9-ca7b-4563-bd44-1918ae6e9f40">

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/f975fd7e-867a-48e9-b63a-5a797af163cf">




