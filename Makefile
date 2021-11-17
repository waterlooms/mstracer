all:
	docker-compose up --build -d

rebuild:
	docker-compose up --build -d --force

shell: force
	docker-compose exec waterlooms bash

force:

