### Building and running your application

First, build your image, e.g.: `docker build -t mt_lintocirc .`.
If your system uses a different CPU architecture than your development
machine (e.g., you are on a Mac M1 and your cloud provider is amd64),
you'll want to build the image for that platform, e.g.:
`docker build --platform=linux/amd64 -t mt_lintocirc .`.

Then, push it to your registry, e.g. `docker push myregistry.com/mt_lintocirc`.

Consult Docker's [getting started](https://docs.docker.com/go/get-started-sharing/)
docs for more detail on building and pushing.

### References
* [Docker's Rust guide](https://docs.docker.com/language/rust/)