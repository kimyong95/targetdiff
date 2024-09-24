aws iam create-role --role-name lambda-ex --assume-role-policy-document '{"Version": "2012-10-17","Statement": [{ "Effect": "Allow", "Principal": {"Service": "lambda.amazonaws.com"}, "Action": "sts:AssumeRole"}]}'

aws lambda create-function --function-name vina-dock --memory-size 1769 --timeout 15 --zip-file fileb://packages.zip --handler lambda_function.lambda_handler --runtime python3.9 --role arn:aws:iam::280819572982:role/lambda-ex

# In AWS console:
# - create function url
# - configure concurrent executions quotas

# rm -r packages/*
cd packages; pip install --platform manylinux2014_x86_64 --target=. --implementation cp --python-version 3.9 --only-binary=:all: vina==1.2.5 numpy==1.23.5 jsonpickle==3.0.4 ; cd ..

rm packages.zip; cd packages; zip -r ../packages.zip * ; cd ..; zip -r packages.zip data/* ; zip -r packages.zip lambda_function.py

# aws lambda update-function-code --function-name vina-dock --zip-file fileb://packages.zip

aws lambda invoke --function-name vina-dock out --log-type Tail --payload file://example_payload.json --cli-binary-format raw-in-base64-out

curl -X POST https://qylhug6f6qwh7xbs5zsewfzwi40hktit.lambda-url.ap-southeast-1.on.aws/ -d @example_payload.json -H "content-type: application/json"