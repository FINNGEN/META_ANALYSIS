pipeline {
   agent any

   stages {
      stage('Build') {
         steps {
            git 'https://github.com/majorseitan/colocalization.git'
            script {  c = docker.build "colocalization" + ":$BUILD_NUMBER"
                      c.inside("-u root"){ sh 'pip install pylint safety pyflakes mypy prospector bandit' }
            }
         }
        }
   }
}
