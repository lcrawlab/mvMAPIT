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
            sshScript remote: remote, script: "install-from-source.sh"
        }
        stage('Submit sbatch job') {
            sshScript remote: remote, script: "sbatch-submit.sh"
        }
    }
}
