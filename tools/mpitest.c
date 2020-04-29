#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int size, rank;

// If MPI's non-overtaking guarantee is met, then the final line should always read "<ABXXXXXXX>", even though the buffered call
// WHOOPS!  I misunderstood the non-overtaking rule!  It's actually more strict (and hence more useful) than I thought:
// I thought that if you were to send 2 messages with different tags and they were received by 2 MPI_ANY_TAG receives, they could
// arrive in either order -- but they must arrive in the order they were sent.  Actually I don't think this fact would have changed
// anything for fastdnamp since at times we need to receive on a *subset* of available tags (i.e. we specifically want to avoid
// receiving certain tags) and MPI provides no way to do this, hence we were required to "pack" several meanings into a single tag
// and receive only on that tag.  But I wonder if we do have any MPI_ANY_TAG receives anywhere...  Could they be screwing things
// up?
void nonovertaking(void) {
	char buf[10] = "XXXXXXXXX";
	MPI_Status statuses[2];
	MPI_Request reqs[2];
	int reqCount;
	int indices[10];
	static unsigned char msgbuf[5000];

	if (rank == 0) {
		strcpy(buf, "ABCDEFG");
		MPI_Buffer_attach(msgbuf, sizeof msgbuf);
		//MPI_Bsend(buf, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		//MPI_Ssend(buf + 1, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		printf("About to send #1.\n");
		MPI_Ibsend(buf, 1, MPI_CHAR, 1, 42, MPI_COMM_WORLD, reqs);
		printf("MPI_Ibsend() #1 returned.\n");
		sleep(2);
		printf("About to send #2.\n");
		MPI_Issend(buf + 1, 1, MPI_CHAR, 1, 69, MPI_COMM_WORLD, reqs + 1);
		printf("MPI_Issend() #2 returned.\n");
		printf("Waiting for sends to complete...\n");
		//MPI_Waitall(2, reqs, statuses);
		MPI_Wait(reqs + 1, statuses + 1);
		MPI_Wait(reqs, statuses);
		MPI_Buffer_detach(msgbuf, (int *) buf);		// Just put the useless results in random locations
	} else {
		sleep(4);
		printf("About to receive #1.\n");
		MPI_Recv(buf, 1, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, statuses);
		printf("About to receive #2.\n");
		MPI_Recv(buf + 1, 1, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, statuses + 1);
		printf("buf=<%s>\n", buf);
	}
}

void pause(int secs) {
	char buf[10];
	MPI_Status status;
	MPI_Request recvReq;
	int reqCount;
	int indices[10];

	if (rank == 0) {
		sleep(secs);
		printf("About to send.\n");
		MPI_Ssend(buf, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		printf("Sent.\n");
	} else {
		printf("About to receive.\n");
		//MPI_Recv(buf, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Irecv(buf, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &recvReq);
		printf("Posted nonblocking receive request.\n");
		//MPI_Wait(&recvReq, &status);
		MPI_Waitsome(1, &recvReq, &reqCount, indices, &status);
		printf("Received.\n");
	}
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (!strcmp(argv[1], "pause")) {
		pause(atoi(argv[2]));
	} else if (!strcmp(argv[1], "nonovertaking")) {
		nonovertaking();
	}

	MPI_Finalize();
	return 0;
}
