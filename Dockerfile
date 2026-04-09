FROM r-base

# Install renv package
RUN R -e 'install.packages("renv", repos="https://cloud.r-project.org")'

RUN mkdir -p /entry
RUN mkdir -p /workspace

WORKDIR /entry
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

RUN chown -R 1000:1000 /entry
RUN chown -R 1000:1000 /workspace
USER 1000

WORKDIR /workspace

ENTRYPOINT ["/entry/entrypoint.sh"]

