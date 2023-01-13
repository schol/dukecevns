from physino/root

run curl -LO https://github.com/schol/dukecevns/archive/refs/heads/master.zip \
 && curl -LO https://github.com/nlohmann/json/archive/refs/heads/develop.zip \
 && unzip master.zip && mv dukecevns-master /root/dukecevns \
 && unzip develop.zip && mv json-develop/* /root/dukecevns/json \
 && rm -fr *.zip json-develop

# overwrite the default working directory: /root/
workdir /root/dukecevns
