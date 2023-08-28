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
<img width="649" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/641216c4-a907-467b-a563-93e8d8fb4339">

<img width="731" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/244a14d4-bd22-433a-b81c-01998f3393cc">

<img width="731" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/f9dadf44-787a-4ee0-98b5-288d9a3fbc1a">

### Value function learning (the state uncertainty + feedback model is shown in the paper)

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/4f8c63f5-c78f-4063-8eca-5520f8a851ee">

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/d20dc79a-a6d5-4d14-ab7c-8aa117e5d1b0">

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/78b1aaf9-ca7b-4563-bd44-1918ae6e9f40">

<img width="526" alt="image" src="https://github.com/harvardschoolofmouse/RPE-Modeling/assets/50119751/f975fd7e-867a-48e9-b63a-5a797af163cf">




