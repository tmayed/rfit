FROM r-base

# 1. Install openssh-server alongside your other dependencies
RUN apt-get update && apt-get install -y curl openssh-server \
    && rm -rf /var/lib/apt/lists/*
    
# Create the privilege separation directory required by SSH
RUN mkdir -p /run/sshd

# 2. Install renv package
RUN R -e 'install.packages("renv", repos="https://cloud.r-project.org")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("dplyr")'
RUN R -e 'install.packages("tidy")'

# 3. Create a dedicated user for ID 1000 (so SSH has a home directory to route to)
# Find the existing user with UID 1000, delete them, and create our appuser
RUN OLD_USER=$(getent passwd 1000 | cut -d: -f1) \
    && if [ -n "$OLD_USER" ]; then userdel -f "$OLD_USER"; fi \
    && useradd -m -u 1000 -s /bin/bash appuser

RUN mkdir -p /entry /workspace

# Copy and setup entrypoint
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

WORKDIR /entry

# Copy and setup entrypoint
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

# Fix permissions
RUN chown -R appuser:appuser /entry
RUN chown -R appuser:appuser /workspace

# IMPORTANT: Do NOT switch to USER 1000 here. 
# We need to stay root so the entrypoint can start the SSH daemon.
# When Gemini SSHs in, it will authenticate as 'appuser' (UID 1000) automatically.

WORKDIR /workspace

ENTRYPOINT ["/entry/entrypoint.sh"]
