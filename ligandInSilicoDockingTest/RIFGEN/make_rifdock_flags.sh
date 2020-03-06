#!/bin/bash

cd $1
tail -n 17 output.log > local_rifdock.flags
cat ./rifdock_v4_ads.flags >> local_rifdock.flags
cd ..
