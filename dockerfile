FROM ubuntu:22.04
WORKDIR /home/anegin97/
# Copy requirements file
COPY requirements.txt /home/anegin97/
RUN ln -fs /usr/share/zoneinfo/Etc/UTC /etc/localtime && \
    echo "Etc/UTC" > /etc/timezone && \
    apt-get update && apt-get install -y tzdata
RUN apt-get update && apt-get install -y r-base-dev
# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    zip \
    unzip \
    wget \
    python3-venv \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libpq-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libssh2-1-dev \
    unixodbc-dev \
    cnvkit \
    libglpk40
# Set up Python virtual environment
RUN python3 -m venv venv
RUN /home/anegin97/venv/bin/pip install --upgrade pip
RUN /home/anegin97/venv/bin/pip install -r /home/anegin97/requirements.txt
# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /home/anegin97/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh
# Initialize Conda for the shell
RUN /home/anegin97/miniconda3/bin/conda init bash
# Ensure Conda doesn't automatically update itself
ENV PATH="/home/anegin97/miniconda3/bin:$PATH"
ENV CONDA_AUTO_UPDATE_CONDA=false
# Create a new Conda environment (if needed)
RUN conda create -n myenv python=3.8 -y
# Complete install guide cnvkit 0.9.10
RUN pip install "numpy>=1.22,<2.1"
RUN wget https://github.com/etal/cnvkit/archive/refs/tags/v0.9.10.zip && \
    unzip v0.9.10.zip && \
    cd cnvkit-0.9.10 && \
    python setup.py install
# Copy all files to the working directory
COPY . /home/anegin97/Bioinformatics/appforclient
# Install Nginx
RUN apt-get update && apt-get install -y nginx && apt-get clean && rm -rf /var/lib/apt/lists/*

# Make shell scripts executable
RUN chmod +x /home/anegin97/Bioinformatics/appforclient/PSeqPGSSE.sh
RUN chmod +x /home/anegin97/Bioinformatics/appforclient/PSeqPGSPE.sh

# Expose port 9999 for the Flask or similar web application
EXPOSE 9999

# Expose port 80 for Nginx (if used)
EXPOSE 80

# Set the command to run Nginx and your Python application
CMD service nginx start && ./venv/bin/python /home/anegin97/Bioinformatics/appforclient/app.py
