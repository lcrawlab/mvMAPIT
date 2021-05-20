#!groovy

node {
    checkout scm
    withCredentials(
            [
                    usernamePassword(
                            credentialsId: 'brown-account',
                            passwordVariable: 'PASSWORD',
                            usernameVariable: 'USER')
            ]
    ) {
        def remote = [
                name: 'oscar',
                host: 'ssh.ccv.brown.edu',
                user: USER,
                password: PASSWORD,
                allowAnyHosts: true
        ]
        stage('Install R package mvMAPIT') {
            sshScript remote: remote, script: "oscar/install-from-source.sh"
        }
        stage('Run R script on login node') {
            sshScript remote: remote, script: "oscar/run-R-script.sh"
        }
        //stage('Submit sbatch job') {
        //    sshScript remote: remote, script: "oscar/sbatch-submit-array.sh"
        //}
    }
}
