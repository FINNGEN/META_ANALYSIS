pipeline {
   agent any

   stages {
      stage('Build') {
         steps {
            script {  c = docker.build "${JOB_NAME}:${BUILD_NUMBER}"
                      c.inside("-u root"){ sh 'pip install pylint safety pyflakes mypy prospector bandit' }
            }
         }
        }
   }
}
