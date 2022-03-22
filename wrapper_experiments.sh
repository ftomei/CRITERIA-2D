num_iterations=150
python src/experiments_launcher.py -nits $num_iterations
python src/json_to_csv.py -nits $num_iterations