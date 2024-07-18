using ReinforcementLearning
run(
    RandomPolicy(),
    CartPoleEnv(),
    StopAfterNSteps(1_000),
    TotalRewardPerEpisode()
)