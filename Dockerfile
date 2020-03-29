FROM continuumio/anaconda3

WORKDIR C:\Users\Fractal\Documents\Github\beanscripts\

COPY . .

# streamlit-specific commands
RUN mkdir -p /root/.streamlit
RUN bash -c 'echo -e "\
[general]\n\
email = \"\"\n\
" > /root/.streamlit/credentials.toml'
RUN bash -c 'echo -e "\
[server]\n\
enableCORS = false\n\
" > /root/.streamlit/config.toml'

# exposing default port for streamlit
EXPOSE 8501

# activate conda environment
SHELL ["/bin/bash", "-c", "source activate ./env"]

CMD [ "streamlit", "run", "C:/Users/Fractal/Documents/Github/beanscripts/first.py" ]