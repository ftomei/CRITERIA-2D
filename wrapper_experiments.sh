conda activate CRITERIA-2D

num_iterations=5
python src/experiments-launcher -nits $num_iterations
python src/json_to_csv -nits $num_iterations